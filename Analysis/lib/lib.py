import uproot
import numpy as np
import pandas as pd
from numba import jit

def combine_data_ophit(data,ophit):
    output = dict()
    output["PE"] = []
    output["TIMES"] = []
    output["AMP"] = []
    for wvf in range(len(data["RECO"]["CH"])):
        ev = data["RECO"]["EV"][wvf]
        ch = data["RECO"]["CH"][wvf]
        ev_ch_filter = (ophit[ev][ch]["PeakTime"] < data["RECO"]["WVF_FX"][wvf]) & (ophit[ev][ch]["PeakTime"] > data["RECO"]["WVF_IX"][wvf])
        output["AMP"].append(ophit[ev][ch]["Amplitude"][ev_ch_filter])
        output["TIMES"].append(ophit[ev][ch]["PeakTime"][ev_ch_filter])
        output["PE"].append(np.sum(ophit[ev][ch]["PE"][ev_ch_filter]))
    return output

def combine_ophit_with_data(data,ophit):
    output = dict()
    output["PE"] = []
    output["TIMES"] = []
    output["AMP"] = []
    for wvf in range(len(data["RECO"]["CH"])):
        ev = data["RECO"]["EV"][wvf]
        ch = data["RECO"]["CH"][wvf]
        output["AMP"].append(ophit[ev][ch]["Amplitude"][np.where((ophit[ev][ch]["PeakTime"] < np.max(data["RECO"]["WVF_X"][wvf])) & (ophit[ev][ch]["PeakTime"] > np.min(data["RECO"]["WVF_X"][wvf])))])
        output["TIMES"].append(ophit[ev][ch]["PeakTime"][np.where((ophit[ev][ch]["PeakTime"] < np.max(data["RECO"]["WVF_X"][wvf])) & (ophit[ev][ch]["PeakTime"] > np.min(data["RECO"]["WVF_X"][wvf])))])
        output["PE"].append(np.sum(ophit[ev][ch]["PE"][np.where((ophit[ev][ch]["PeakTime"] < np.max(data["RECO"]["WVF_X"][wvf])) & (ophit[ev][ch]["PeakTime"] > np.min(data["RECO"]["WVF_X"][wvf])))]))
    return output

def ophit_data(wvf_path,scaling,label,branches=[],as_df=False,debug=False):
    ophit = dict()
    ophit["SCALING"] = scaling
    ophit["LABEL"] = label
    with uproot.open(wvf_path) as raw:
        print(raw["ophitana/PerOpHitTree;1"].keys())
        if branches == []: branch_list = raw["ophitana/PerOpHitTree;1"].keys()
        else: branch_list = branches
        type_dict = {"EventID":np.int8,"OpChannel":np.int16,"PeakTime":np.float32,"Amplitude":np.float32,"PE":np.float32}
        for branch in branches:
            ophit[branch] = raw["ophitana/PerOpHitTree;1"][branch].array().to_numpy().astype(type_dict[branch])
    if as_df == True:
        return pd.DataFrame(ophit,columns=branch_list)
    else:
        return ophit

def order_ophit_data(ophit,max_ev,ch_array):    
    new_ophit = dict()
    new_ophit["SCALING"] = ophit["SCALING"]
    new_ophit["LABEL"] = ophit["LABEL"]
    for ev in range(max_ev):
        new_ophit[ev+1] = dict()
        for ch in range(np.max(ch_array)+1):
            new_ophit[ev+1][ch] = dict()
            for key in ophit.keys():
                if key != "LABEL" and key != "SCALING":
                    new_ophit[ev+1][ch][key] = ophit[key][(ophit["OpChannel"] == ch)*(ophit["EventID"] == ev+1)]
    return new_ophit

def get_wvf_data(wvf_type,wvf_path,root_folder,max_ev,pedestal,output,template,debug=False):
    with uproot.open(wvf_path) as raw:
        raw_wvf_list = raw[root_folder[1]].keys()
        wvf_idx_map = {h: i for i, h in enumerate(raw_wvf_list)}
        filtered_hist = [h for h in raw_wvf_list if "event" in h and int(h.split("_")[1]) <= max_ev]

        idx_list = [wvf_idx_map[h] for h in filtered_hist]
        if wvf_type == "RAW":
            pe_per_wvf = get_true_pe(raw, root_folder, max_ev, debug=debug)
            print(pe_per_wvf.keys())
            raw_wvfs_pe = np.empty(len(idx_list), dtype=object)
            raw_wvfs_pe_count = np.empty(len(idx_list), dtype=np.int32)
            raw_wvfs_t0 = np.empty(len(idx_list), dtype=np.float32)

        raw_wvfs =       np.empty([len(idx_list),1000], dtype=np.float32)
        raw_wvfs_ix =    np.empty(len(idx_list), dtype=np.float32)
        raw_wvfs_fx =    np.empty(len(idx_list), dtype=np.float32)
        raw_event =      np.empty(len(idx_list), dtype=np.int8)
        raw_wvf_ch =     np.empty(len(idx_list), dtype=np.int16)
        raw_wvf_ch_num = np.empty(len(idx_list), dtype=np.int8)

        for idx,raw_idx in enumerate(idx_list):
            this_raw_wvf   = raw[root_folder[1]][raw_wvf_list[raw_idx]].to_numpy()[0]
            this_raw_wvf_x = raw[root_folder[1]][raw_wvf_list[raw_idx]].to_numpy()[1]
            
            if this_raw_wvf.size <= 1000:
                raw_wvfs[idx] = np.pad(this_raw_wvf, (0,1000 - this_raw_wvf.size), 'constant', constant_values=(output["PEDESTAL"],output["PEDESTAL"]))
            else:
                raw_wvfs[idx] = this_raw_wvf[:1000]

            raw_wvfs_ix[idx] = this_raw_wvf_x[0]+float(raw_wvf_list[raw_idx].split("_")[7])
            raw_wvfs_fx[idx] = this_raw_wvf_x[-1]+float(raw_wvf_list[raw_idx].split("_")[7])
            raw_event[idx] = int(raw_wvf_list[raw_idx].split("_")[1])
            raw_wvf_ch[idx] = int(raw_wvf_list[raw_idx].split("_")[3])
            raw_wvf_ch_num[idx] = int(raw_wvf_list[raw_idx].split("_")[5])
            
            if wvf_type == "RAW":
                this_ch_pe = pe_per_wvf.loc[(pe_per_wvf["EV"] == raw_event[idx])*(pe_per_wvf["CH"] == raw_wvf_ch[idx])]["PE"].to_numpy()
                raw_wvfs_pe[idx] = this_ch_pe[(this_ch_pe > raw_wvfs_ix[idx]) & (this_ch_pe < raw_wvfs_fx[idx])]
                raw_wvfs_pe_count[idx] = raw_wvfs_pe[idx].size
                try:
                    raw_wvfs_t0[idx] = raw_wvfs_pe[idx][0]
                except:
                    raw_wvfs_t0[idx] = 0

    output["RECO"]["EV"]     = raw_event    
    output["RECO"]["CH"]     = raw_wvf_ch                 
    output["RECO"]["#WVF"]   = raw_wvf_ch_num                 
    output["RECO"]["WVF"]    = raw_wvfs
    output["RECO"]["WVF_IX"] = raw_wvfs_ix
    output["RECO"]["WVF_FX"] = raw_wvfs_fx
    output["RECO"]["T0"]     = output["SAMPLING"]*1e6*np.argmax(raw_wvfs,axis=1)+raw_wvfs_ix
    output["RECO"]["AMP"]    = np.max(np.asarray(raw_wvfs),axis=1) - output["PEDESTAL"]

    if wvf_type == "RAW": 
        output["RECO"]["PE"]     = np.sum(raw_wvfs-output["PEDESTAL"],where=raw_wvfs-output["PEDESTAL"] > 0,axis=1)/np.sum(template[template > 0])      
        output["TRUE"]["PETIMES"] = raw_wvfs_pe
        output["TRUE"]["PE"] = raw_wvfs_pe_count
        output["TRUE"]["T0"] = raw_wvfs_t0
    if wvf_type == "DEC":
        output["RECO"]["PE"] = np.sum(raw_wvfs-output["PEDESTAL"],where=raw_wvfs-output["PEDESTAL"] > 0,axis=1)
    return output

def get_true_pe(raw,root_folder,max_ev,debug=False):
    # Include true photon information
    true_ev = np.empty(1, dtype=np.int8)  
    true_ch = np.empty(1, dtype=np.int16)  
    true_pe = np.empty(1, dtype=np.float32)

    for ev_idx in range(max_ev):
        try:
            this_ev = raw[root_folder[0]]["PhotonData;1"]["EventID"].array()[ev_idx]
        except:
            this_ev = int(ev_idx+1)

        if this_ev != ev_idx+1:
            print("ERROR: Event ID mismatch")
            print("Event ID: {}".format(this_ev))
            print("Event ID: {}".format(ev_idx+1))
        
        true_ev = np.concatenate((true_ev,this_ev*np.ones(len(raw[root_folder[0]]["PhotonData;1"]["photon_opCh"].array()[ev_idx]))))
        true_ch = np.concatenate((true_ch,raw[root_folder[0]]["PhotonData;1"]["photon_opCh"].array()[ev_idx]))
        true_pe = np.concatenate((true_pe,raw[root_folder[0]]["PhotonData;1"]["photon_pulse"].array()[ev_idx]))
        pe_per_wvf = pd.DataFrame(data={"EV":true_ev,"CH":true_ch},dtype=np.int16)
        pe_per_wvf["PE"] = true_pe
        
    return pe_per_wvf    

def load_root_new(wvf_type,wvf_path,root_folder,template,label="",max_ev=1,debug=False):
    # Generate empty output dict
    output = dict()
    output["SAMPLING"] = 16e-9
    output["RECO"] = dict()
    # Load true photon data
    if wvf_type == "RAW":
        output["TRUE"] = dict()
        output["PE_INFO"] = dict()
        output["PEDESTAL"] = 1500
        if label == "": output["LABEL"] = "RAW"
        else: output["LABEL"] = label

    # Set pedestal value
    else: output["PEDESTAL"] = 0
    # Load waveform data
    get_wvf_data(wvf_type,wvf_path,root_folder,max_ev,output["PEDESTAL"],output,template,debug=debug)
    if wvf_type == "DEC":
        if label == "": output["LABEL"] = "DEC"
        else: output["LABEL"] = label
    return output
    
def load_larsoft_root(wvf_type,wvf_path,root_folder,template,label="",max_ev=1,debug=False):
    # Generate empty output dict
    output = dict()
    output["RECO"] = dict()

    # Open _hist.root data
    with uproot.open(wvf_path) as raw:
        raw_wvf_list = raw[root_folder[1]].keys()
        raw_wvfs           = []
        raw_wvfs_x         = []
        raw_event          = []
        raw_wvf_ch         = []
        raw_wvf_ch_num     = []
        raw_wvf_pe         = []
        raw_wvf_first_time = []
        raw_wvf_time       = []
        
        # Loop over histograms in wvf data folder
        for this_wvf in raw_wvf_list:
            # Avoid extra objects in hist folder
            if this_wvf == "TriggerData;1" or this_wvf == "Daphne;1":
                continue
            else:
                # Extract information from hist title
                wvf_info = this_wvf.split("_")
                wvf_info[-1] = wvf_info[-1].split(";")[0]
                
                if int(wvf_info[1]) <= max_ev:
                    raw_event.append(int(wvf_info[1]))
                    raw_wvf_ch.append(int(wvf_info[3]))
                    raw_wvf_ch_num.append(int(wvf_info[5]))
                    raw_wvf_first_time.append(float(wvf_info[7]))
                    raw_wvf_time.append(float(wvf_info[9]))
                    
                    # Save wvf data (y and x axis)
                    raw_wvfs.append(raw[root_folder[1]][this_wvf].to_numpy()[0])
                    raw_wvfs_x.append(float(wvf_info[7])+raw[root_folder[1]][this_wvf].to_numpy()[1][:-1])
                else:
                    continue

        output["RECO"]["EV"]    = np.asarray(raw_event)    
        output["RECO"]["CH"]    = np.asarray(raw_wvf_ch)                 
        output["RECO"]["#WVF"]  = np.asarray(raw_wvf_ch_num)                 
        output["RECO"]["WVF"]   = raw_wvfs
        output["RECO"]["WVF_X"] = raw_wvfs_x

        for key in ["PE","AMP","PETIMES","T0","PED","PEDSTD"]:
            output["RECO"][key] = []
        
        if wvf_type == "RAW":
            output["TRUE"] = dict()
            output["PRETRIGGER"] = 100
            output["PEDESTAL"] = 1500
            if label == "": output["LABEL"] = "RAW"
            else: output["LABEL"] = label
            
            # Include true photon information
            output["TRUTH_EV"] = np.empty(1)  
            output["TRUTH_CH"] = np.empty(1)  
            output["TRUTH_PE"] = np.empty(1)

            for ev_idx in range(max_ev):
                try:
                    this_ev = raw[root_folder[0]]["PhotonData;1"]["EventID"].array()[ev_idx]
                except:
                    this_ev = ev_idx+1

                if this_ev != ev_idx+1:
                    print("ERROR: Event ID mismatch")
                    print("Event ID: {}".format(this_ev))
                    print("Event ID: {}".format(ev_idx+1))

                output["TRUTH_EV"] = np.concatenate((output["TRUTH_EV"],this_ev*np.ones(len(raw[root_folder[0]]["PhotonData;1"]["photon_opCh"].array()[ev_idx]))))
                output["TRUTH_CH"] = np.concatenate((output["TRUTH_CH"],raw[root_folder[0]]["PhotonData;1"]["photon_opCh"].array()[ev_idx]))
                output["TRUTH_PE"] = np.concatenate((output["TRUTH_PE"],raw[root_folder[0]]["PhotonData;1"]["photon_pulse"].array()[ev_idx]))

            output["TRUE"]["PETIMES"]  = []
            output["TRUE"]["PE"]  = []
            output["TRUE"]["T0"]  = []
            
            photons_per_wvf = photon_arrival_times(output,max_ev)

            for i in range(len(output["RECO"]["WVF"])):
                this_ev = output["RECO"]["EV"][i]
                
                try:output["TRUE"]["PE"].append(len(photons_per_wvf[this_ev][i]))
                except:output["TRUE"]["PE"].append(0)

                try:output["RECO"]["PE"].append(np.sum(np.asarray(raw_wvfs[i]-output["PEDESTAL"])[raw_wvfs[i]-output["PEDESTAL"] > 0])/np.sum(template[template > 0]))
                except:output["RECO"]["PE"].append(0)
                
                try:output["TRUE"]["PETIMES"].append(photons_per_wvf[this_ev][i])
                except:output["TRUE"]["PETIMES"].append([])

                try:output["TRUE"]["T0"].append(photons_per_wvf[this_ev][i][0])
                except IndexError:output["TRUE"]["T0"].append(0)
                except ValueError:output["TRUE"]["T0"].append(0)

        elif wvf_type == "DEC":
            output["PRETRIGGER"] = 100
            output["PEDESTAL"] = 0
            if label == "": output["LABEL"] = "DEC"
            else: output["LABEL"] = label
            for i in range(len(output["RECO"]["WVF"])):
                output["RECO"]["PE"].append(np.sum(np.asarray(raw_wvfs[i]-output["PEDESTAL"])))

        else:
            output["PRETRIGGER"] = 0
            output["PEDESTAL"] = 0

        for i in range(len(output["RECO"]["WVF"])):
            output["RECO"]["T0"].append(raw_wvfs_x[i][np.argmax(raw_wvfs[i])])
            output["RECO"]["AMP"].append(np.max(raw_wvfs[i])-output["PEDESTAL"])      
            output["RECO"]["PED"].append(np.mean(raw_wvfs[i][:output["PRETRIGGER"]])-output["PEDESTAL"])      
            output["RECO"]["PEDSTD"].append(np.std(raw_wvfs[i][:output["PRETRIGGER"]])-output["PEDESTAL"])

        output["SAMPLING"] = 16e-9

    return output

def photon_arrival_times(raw,max_ev):
    photons_per_wvf = dict()
    if np.max(raw["RECO"]["EV"]) <= max_ev:
        max_ev = np.max(raw["RECO"]["EV"])
    
    for jj in range(len(raw["RECO"]["WVF_X"])):
        ev = raw["RECO"]["EV"][jj]
        if ev > max_ev:
            break

        min_wvf_time = np.min(raw["RECO"]["WVF_X"][jj])
        max_wvf_time = np.max(raw["RECO"]["WVF_X"][jj])
        photons_per_channel = raw["TRUTH_PE"][(raw["TRUTH_CH"] == raw["RECO"]["CH"][jj])*(raw["TRUTH_EV"] == ev)]
        try:
            photons_per_wvf[ev].append(photons_per_channel[(photons_per_channel > min_wvf_time) & (photons_per_channel < max_wvf_time)])
        except KeyError:
            photons_per_wvf[ev] = [photons_per_channel[(photons_per_channel > min_wvf_time) & (photons_per_channel < max_wvf_time)]]

    return photons_per_wvf

def order_wvfs(wvfs,key_list):
    output = dict()
    for ii in range(np.max(wvfs["RECO"]["EV"])):
        output[ii+1] = dict()
    
    for jj in range(len(wvfs["RECO"]["CH"])):
        output[wvfs["RECO"]["EV"][jj]][wvfs["RECO"]["CH"][jj]] = dict()
        for key in key_list:
            output[wvfs["RECO"]["EV"][jj]][wvfs["RECO"]["CH"][jj]][key] = []

    for kk in range(len(wvfs["RECO"]["WVF"])):
        for key in key_list:
            output[wvfs["RECO"]["EV"][kk]][wvfs["RECO"]["CH"][kk]][key].append(wvfs["RECO"][key][kk])
    
    return output

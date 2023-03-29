import uproot
import numpy as np

def scint_profile(x,const,a_f,tau_f,tau_s):
    return const*(2*a_f/tau_f*np.exp(-(x)/tau_f) + 2*(1-a_f)/tau_s*np.exp(-(x)/tau_s))

def combine_ophit_with_data(data,ophit):
    output = dict()
    output["PE"] = []
    output["TIMES"] = []
    output["AMP"] = []
    for wvf in range(len(data["RECO"]["CH"])):
        ch = data["RECO"]["CH"][wvf]
        output["AMP"].append(ophit[ch]["Amplitude"][np.where((ophit[ch]["PeakTime"] < np.max(data["RECO"]["WVF_X"][wvf])) & (ophit[ch]["PeakTime"] > np.min(data["RECO"]["WVF_X"][wvf])))])
        output["TIMES"].append(ophit[ch]["PeakTime"][np.where((ophit[ch]["PeakTime"] < np.max(data["RECO"]["WVF_X"][wvf])) & (ophit[ch]["PeakTime"] > np.min(data["RECO"]["WVF_X"][wvf])))])
        output["PE"].append(np.sum(ophit[ch]["PE"][np.where((ophit[ch]["PeakTime"] < np.max(data["RECO"]["WVF_X"][wvf])) & (ophit[ch]["PeakTime"] > np.min(data["RECO"]["WVF_X"][wvf])))]))
    return output

def ophit_data(wvf_path,scaling,label):
    ophit = dict()
    ophit["SCALING"] = scaling
    ophit["LABEL"] = label
    raw = uproot.open(wvf_path)
    print(raw["ophitana/PerOpHitTree;1"].keys())
    for branch in raw["ophitana/PerOpHitTree;1"].keys():
        ophit[branch] = raw["ophitana/PerOpHitTree;1"][branch].array()
    return ophit

def order_ophit_data(ophit,ch_array):
    new_ophit = dict()
    new_ophit["SCALING"] = ophit["SCALING"]
    new_ophit["LABEL"] = ophit["LABEL"]
    for ch in ch_array:
        new_ophit[ch] = dict()
        for key in ophit.keys():
            if key != "LABEL" and key != "SCALING":
                new_ophit[ch][key] = ophit[key][np.where(ophit["OpChannel"] == ch)]
    return new_ophit

def load_larsoft_root(wvf_type,wvf_path,root_folder,template,label="",debug=False):
    # Generate empty output dict
    output = dict()
    output["RECO"] = dict()

    # Open _hist.root data
    raw = uproot.open(wvf_path)
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
            raw_event.append(int(wvf_info[1]))
            raw_wvf_ch.append(int(wvf_info[3]))
            raw_wvf_ch_num.append(int(wvf_info[5]))
            raw_wvf_first_time.append(float(wvf_info[7]))
            raw_wvf_time.append(float(wvf_info[9]))
            
            # Save wvf data (y and x axis)
            raw_wvfs.append(raw[root_folder[1]][this_wvf].to_numpy()[0])
            raw_wvfs_x.append(float(wvf_info[7])+raw[root_folder[1]][this_wvf].to_numpy()[1][:-1])

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
        raw_ch = raw[root_folder[0]]["PhotonData"]["photon_opCh"].array().to_numpy()[0]
        raw_npe = raw[root_folder[0]]["PhotonData"]["photon_pulse"].array().to_numpy()[0]
        output["TRUTH_CH"]  = raw_ch   
        output["TRUTH_PE"]  = raw_npe
        output["TRUE"]["PETIMES"]  = []
        output["TRUE"]["PE"]  = []
        output["TRUE"]["T0"]  = []

        photons_per_wvf = photon_arrival_times(output)
        for i in range(len(output["RECO"]["WVF"])):
            output["TRUE"]["PE"].append(len(photons_per_wvf[i]))
            output["RECO"]["PE"].append(np.sum(np.asarray(raw_wvfs[i]-output["PEDESTAL"])[raw_wvfs[i]-output["PEDESTAL"] > 0])/np.sum(template[template > 0]))

            # try:
            output["TRUE"]["PETIMES"].append(photons_per_wvf[i])
            try:
                output["TRUE"]["T0"].append(photons_per_wvf[i][0])
            except IndexError:
                output["TRUE"]["T0"].append(0)

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

def photon_arrival_times(raw):
    photons_per_wvf = []

    for jj in range(len(raw["RECO"]["CH"])):
        min_wvf_time = np.min(raw["RECO"]["WVF_X"][jj])
        max_wvf_time = np.max(raw["RECO"]["WVF_X"][jj])
        photons_per_channel = raw["TRUTH_PE"][np.where(raw["TRUTH_CH"] == raw["RECO"]["CH"][jj])]
        photons_per_wvf.append(photons_per_channel[(photons_per_channel > min_wvf_time) & (photons_per_channel < max_wvf_time)])
    
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

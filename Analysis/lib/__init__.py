import os

import numpy as np
import matplotlib.pyplot as plt
import plotly.figure_factory as ff
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd

from scipy.optimize  import curve_fit
from plotly.subplots import make_subplots
from itertools       import product

from .lib  import *
from .fit  import *
from .plot import *

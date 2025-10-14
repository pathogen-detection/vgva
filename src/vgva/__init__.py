import datetime
import matplotlib as mpl
import Bio
from .simi.simi import Similarity
from .genv.genome_variation import GenomeVariation
from .stat.statistical_analysis import calculate_boxplot_parameters

assert Bio.__version__ == '1.83', "biopython requires BioPython v1.8 or newer"
__author__ = "Bangjun Gong"
__email__ = "bangjungong@foxmail.com"
__version__ = "1.0.0"
__license__ = "GPL v3"
__copyright__ = f"Copyright (c) {datetime.datetime.now().year} Bangjun Gong"


def set_text_editable(_bool=True):
    if not _bool:
        return
    mpl.rcParams['svg.fonttype'] = 'none'

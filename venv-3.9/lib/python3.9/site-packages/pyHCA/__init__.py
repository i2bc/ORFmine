__version__ = "0.5"
__release__ = __version__  + '-dev' # comment out '-dev' before a release


import sys
from .core.annotateHCA import main_segment as annotate_main
from .core.domainseq_segHCA import main_domainseq as domainseq_main
from .core.drawHCA import main as draw_main
from .core.tremoloHCA  import main as tremolo_main
from .core.HCA import HCA
from .core.disorderHCA import main as disorder_main

ete3_available = True
dom_on_tree_available = True
try:
    from .core.domains_on_tree  import main as dom_ontree_main
except:
    dom_on_tree_available = False
    try:
        from ete3 import SeqMotifFace
    except ImportError:
        ete3_available = False
        print('An error occured while importing ete3, please check that ete3 is correctly installed allong PyQt4 (run python -c "from ete3 import SeqMotifFace"', file=sys.stderr)
    print("Something went wrong importing dom_ontree_main", file=sys.stderr)

#print(dom_on_tree_available)


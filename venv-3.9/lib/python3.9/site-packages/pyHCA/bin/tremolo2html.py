#!/usr/bin/env python
""" script to transform tremolo text result file into an html version for easier visualization
"""

import os, sys, argparse
import numpy as np
from pyHCA.core.drawHCA import createHCAsvg
from pyHCA.core.ioHCA import read_tremolo
from pyHCA.core.annotateHCA import _annotation_aminoacids as segmentation
from myPython.Residu import AA1

def get_cmd():
    """ get command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", action="store", dest="inputres", help="input tremolo result file", required=True)
    parser.add_argument("-q", action="store", dest="query", help="query sequence", required=True)
    parser.add_argument("-w", action="store", dest="workdir", help="working directory to store html data", required=True)
    parser.add_argument("-o", action="store", dest="output", help="html output file", required=True)
    parser.add_argument("--ffdata", action="store", dest="ffdata", help="fasta file", required=True)
    parser.add_argument("--ffindex", action="store", dest="ffindex", help="ffindex file", required=True)
    params = parser.parse_args()
    return params

### HTML 

def html_header():
    """
    brief add DOCTYPE and header and script 
    """
    
    return """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
    
    <header>
        <style type="text/css"><!--
            .alig { font-family: "Courier New", Courier; font-size:8pt; }
            html, body { margin:0; padding:0 }
            object { width: 100%; height: 100% }
            --></style>\n\n
                    
                    
        <script language="javascript">
            function showhide(id, svg) {
                var element = document.getElementById(id);
                if (element.style.display == "none") {
                    // Check if SVG object already loaded; if not, load it now
                    if (!element.getElementsByTagName("object").length) {
                        var object = document.createElement("object");
                        object.type = "image/svg+xml";
                        object.data = svg;
                        element.appendChild(object);
                    }
                    element.style.display = "block";
                }
                else {
                    element.style.display = "none";
                }
            }
            function loadSVG(id, svg) {
                var element = document.getElementById(id);
                var object = document.createElement("object");
                object.type = "image/svg+xml";
                object.data = svg;
                element.appendChild(object);
            }
                
        </script>
    </header>
    <body>
"""

def html_footer():
    """ footer 
    """
    return '</body>\n</html>\n'

def create_protheader(prot):
    """ write protein header
    """
    name = prot.split("|")[1]
    return '        <li><a href="http://www.uniprot.org/uniprot/{}">>{}</a>\n'.format(prot.split("|")[1], prot)

def create_hitinfo(data):
    """ create hit info
    """
    return "            <li><b>Eval={} ; Probab={} ; Id={} ; Sim={} ; </b></br>\n".format(data["E-value"], data["Probab"], data["Identities"], data["Similarity"])

### SVG
def svg_tremolo(dfasta, hits, domains, dhca, sizes, workdir):
    """ draw tremolo results
    """
    pathsvg = os.path.join(workdir, "svg")
    if not os.path.isdir(pathsvg):
        os.makedirs(pathsvg)
        
    writeAjax(pathsvg)
    writeCSS(pathsvg)
    
    # for each protein draw its domain arrangement, its hydrophobic cluster domains, its clusters, its query hits
    # also draw the alignment between 

    for prot in dhca:
        if prot != "query":
            create_annotation(prot, hits[prot], sizes[prot], domains[prot], dhca[prot], pathsvg)
    
    # create alignment
    create_alitopo(hits, dhca, dfasta, pathsvg)
    
    # create hca svg
    create_hcasvg(hits, dfasta, pathsvg)

# SVG ANNOTATION
def svg_hca_header(nbaa):
    return """<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg width="{}" height="100" version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">

 <!-- ECMAScript with each click -->
  <script type="application/ecmascript"> <![CDATA[
    function aa_click(evt) {{
      var aa = evt.target;
      var currentColor = aa.getAttribute("style");
      if (currentColor == "fill:blue; fill-opacity:0.0")
        {{aa.setAttribute("style", "fill:blue; fill-opacity:0.3");
        }}
      if (currentColor == "fill:blue; fill-opacity:0.3")
        {{aa.setAttribute("style", "fill:red; fill-opacity:0.3");
        }}
      if (currentColor == "fill:red; fill-opacity:0.3")
        {{aa.setAttribute("style", "fill:gray; fill-opacity:0.3");
        }}
      if (currentColor == "fill:gray; fill-opacity:0.3")
        {{aa.setAttribute("style", "fill:blue; fill-opacity:0.0");
        }}
        
    }}
  ]]> </script>
""".format(((nbaa)/20.0)*90+150)

def create_hcasvg(dhits, dfasta, pathsvg):
    """ create hca plot for hits and query sub sequences
    """
    qseq = dfasta["query"]
    for prot in dhits:
        if prot not in dfasta:
            continue
        tseq = dfasta[prot]
        for hitnum in dhits[prot]:
            # get start, stop and sub sequences
            qpathout = os.path.join(pathsvg, "hca_q_{}_{}.svg".format(prot, hitnum)) 
            tpathout = os.path.join(pathsvg, "hca_t_{}_{}.svg".format(prot, hitnum))
            qstart, qstop = dhits[prot][hitnum]["Qstart"], dhits[prot][hitnum]["Qstop"]
            tstart, tstop = dhits[prot][hitnum]["Tstart"], dhits[prot][hitnum]["Tstop"]
            qsubseq = qseq[qstart: qstop]
            tsubseq = tseq[tstart: tstop]
            
            with open(tpathout, "w") as outf:
                tsvg = createHCAsvg(tsubseq, len(tsubseq), domains=[], b=0, yoffset=0)
                # header
                outf.write(svg_hca_header(len(tsubseq)))
                outf.write(tsvg)
                outf.write("</svg>")
            with open(qpathout, "w") as outf:
                qsvg = createHCAsvg(qsubseq, len(qsubseq), domains=[], b=0, yoffset=0)
                # header
                outf.write(svg_hca_header(len(qsubseq)))
                outf.write(qsvg)
                outf.write("</svg>")
    
                
def create_domain(domain, start, stop):
    """ create domain
    """
    color = "gray"
    name = "-".join([domain, str(start), str(stop)])
    return '<rect id="{}" x="{}" y="{}" width="{}" height="10" onmousemove="ShowTooltip(evt, \'{}\')" onmouseout="HideTooltip(evt)" style="fill:{};stroke-width:2;stroke:black" />\n'.format(name, start, 17.5, stop-start+1, domain, color)

def create_annotation(prot, hits, size, domains, dhca, workdir):
    """ create svg figure corresponding to domain _annotation_aminoacids
    """
    pathout = os.path.join(workdir, prot+".svg")
    
    with open(pathout, "w") as outf:
        # svg header
        outf.write("""<?xml version="1.0" encoding="UTF-8"?>
    <?xml-stylesheet href="style.css" type="text/css"?>
    <!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">""")

        outf.write('<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="{}" height="50" onload="init(evt)">\n\n'.format(size+50))
        outf.write('<script xlink:href="script.js" type="text/ecmascript"/>')
        
        # protein sequence
        outf.write('<rect x="0" y="20"   width="{}" height="5"/>\n\n'.format(size))
    
        # protein amino acid scale
        for j in range(50, int(size), 50):
            outf.write('<text x="{}" y="10" style="font-size:7px;">{}</text>\n'.format(j, j))
        
        # draw each domains
        for start, stop, domain in domains:
            line = create_domain(domain, start, stop)
            outf.write(line)
        
        # hca domains
        for dom in dhca["domain"]:
            begin, end = dom.start, dom.stop
            name = "HHCD.{}-{}".format(begin, end)
            outf.write('<rect id="{}" x="{}" y="{}" width="{}" height="10" onmousemove="ShowTooltip(evt, \'{}\')" onmouseout="HideTooltip(evt)" style="fill:grey;" />\n'.format(name, begin, 30, end-begin+1, name))
            
        # hca amas
        for clust in dhca["cluster"]:
            begin, end = clust.start, clust.stop
            outf.write('<rect id="#" x="{}" y="{}" width="{}" height="3" style="fill:black;stroke-width:0" />\n'.format(begin, 30, end-begin+1,))
            
        # hits
        ypos = 0
        for hitnum in hits:
            begin, end = hits[hitnum]["Tstart"], hits[hitnum]["Tstop"]
            evalue = hits[hitnum]["E-value"]
            
            if float(evalue) < 0.005:
                color = "red"
            else:
                color = "orange"
                
            annot = "{}-{} {}".format(begin, end, evalue)
            outf.write('<line x1="{}" y1="{}" x2="{}" y2="{}" onmousemove="ShowTooltip(evt, \'{}\')" onmouseout="HideTooltip(evt)" style="fill:{};stroke:red;stroke-width:3;stroke-linecap:round;" />'.format(begin, ypos+10, end, ypos+10, annot, color))
            
            # if too many hit we up the hit in the svg 
            if ypos < 50:
                ypos += 2
            else:
                ypos -= 2
        
        # footer of svg elements
        outf.write(svg_footer())
    
def svg_footer():
    """ add hidden panel to work with ECMAscript
    """
    
    return """<rect class="tooltip_bg" id="tooltip_bg"
x="0" y="0" rx="4" ry="4"
width="55" height="17" visibility="hidden"/>
<text class="tooltip" id="tooltip"
x="0" y="0" visibility="hidden">Tooltip</text>
</svg>"""

# SVG ALIGNMENT

    
def create_alitopo(dhits, dhca, dfasta, pathsvg):
    """ create svg file containing alignment of query and hit
    """
    letter_spacing = 7
    line_spacing = 50
    letter_begin = 80
    
    yquery = 20
    subject_spacing = 35
    
    masteramas = set()
    for clust in dhca["query"]["cluster"]:
        masteramas.update(set(range(clust.start, clust.stop)))
        
    for prot in dhits:
        amas = set()
        if prot not in dhca:
            continue
        for clust in dhca[prot]["cluster"]:
            amas.update(set(range(clust.start, clust.stop)))
            
        for hitnum in dhits[prot]:
            pathout = os.path.join(pathsvg, "{}_{}.svg".format(prot, hitnum))
            qstart, qstop = dhits[prot][hitnum]["Qstart"], dhits[prot][hitnum]["Qstop"]
            tstart, tstop = dhits[prot][hitnum]["Tstart"], dhits[prot][hitnum]["Tstop"]
            qseqali, tseqali = dhits[prot][hitnum]["Qali"], dhits[prot][hitnum]["Tali"]
            # list of position in the hit of the query
            listposition = range(qstart, qstop+1)
            itepos = 0
            nbconservedtopo = 0
            hydrophobe = "YIMLFWV"
            
            with open(pathout, "w") as outf:
                # header
                
                outf.write('<?xml version="1.0" standalone="no"?>\n')
                outf.write('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">')
                outf.write('<svg xmlns="http://www.w3.org/2000/svg" width="800" height="{}" font-family="Courier New">\n\n'.format((((qstop-qstart)//80)+1)*30))
                for iteline, offset in enumerate(range(0, len(qseqali), 80)):
                    qbegin, qend = qstart + 80 * iteline, min(qstart + 80 + (80 * iteline), qstop)
                    tbegin, tend = tstart + 80 * iteline, min(tstart + 80 + (80 * iteline), tstop)
                    posstop = 700
                    if len(qseqali) < 80:
                        posstop = 7.5 * len(qseqali) + letter_begin

                    outf.write('<text x="10" y="{}" style="font-size:12px;">{:6}</text>\n'.format(yquery+iteline*line_spacing, qbegin))
                    outf.write('<text x="{}" y="{}" style="font-size:12px;">{:6}</text>\n'.format(posstop, yquery+iteline*line_spacing, qstop))
                    outf.write('<text x="10" y="{}" style="font-size:12px;">{:6}</text>\n'.format(subject_spacing+iteline*line_spacing, tbegin))
                    outf.write('<text x="{}" y="{}" style="font-size:12px;">{:6}</text>\n'.format(posstop, subject_spacing+iteline*line_spacing, tend))
                    
                    qali = qseqali[offset: offset+80]
                    tali = tseqali[offset: offset+80]
                    # aa in the line
                    dite2subject = {}
                    current_subject_position = offset
                    for ite, aa in enumerate(tali):
                        if aa != "-":
                            dite2subject[ite] = current_subject_position
                            current_subject_position += 1 
                        else:
                            dite2subject[ite] = ""
                    for itequery, aa in enumerate(qali):
                        if aa != "-":
                            # get the real position in the query sequence
                            cur_position = listposition[itepos]
                            itepos += 1
                            #if cur_position not in dpos2topo:
                                #print("ERROR no match at this position!!!!!")
                                #sys.exit("ERROR colorAlig function")
                            # key are isqueryhydro ishydro istopo
                            #topohydro = dpos2topo[cur_position]
                            #istopo = topohydro["istopo"]
                            #ishydro = topohydro["ishydro"]
                            istopo = False
                            ishydro = False
                        else:
                            istopo = False
                            ishydro = False
                        #if dite2subject[itequery] in amas:
                            #outf.write('<rect x="{}" y="{}" width="7" height="2" style="fill:black; fill-opacity:1.0;" />\n'.format(letter_begin+itequery*letter_spacing, subject_spacing+iteline*line_spacing-12))
                        #if cur_position in masteramas:
                            #outf.write('<rect x="{}" y="{}" width="7" height="2" style="fill:black; fill-opacity:1.0;" />\n'.format(letter_begin+itequery*letter_spacing, yquery+iteline*line_spacing-12))
                        # <Hydrophobe Alignment>
                        if ishydro:
                            # <Topohydrophobe>
                            if istopo:
                                # <SVG::query>
                                outf.write('<rect x="{}" y="{}" width="7" height="10" style="fill:blue; fill-opacity:0.5;" />\n'.format(letter_begin+itequery*letter_spacing, yquery+iteline*line_spacing-10))
                                outf.write('<text x="{}" y="{}" style="font-size:12px;">{}</text>\n'.format(letter_begin+itequery*letter_spacing, yquery+iteline*line_spacing, aa))
                                if tali[itequery] in hydrophobe:
                                    nbconservedtopo += 1
                                    # <SVG::subject>
                                    outf.write('<rect x="{}" y="{}" width="7" height="10" style="fill:blue; fill-opacity:0.5;" />\n'.format(letter_begin+itequery*letter_spacing, subject_spacing+iteline*line_spacing-10))
                                    outf.write('<text x="{}" y="{}" style="font-size:12px;">{}</text>\n'.format(letter_begin+itequery*letter_spacing, subject_spacing+iteline*line_spacing, tali[itequery]))
                                else:
                                    # <SVG::subject>
                                    outf.write('<rect x="{}" y="{}" width="7" height="10" style="fill:red; fill-opacity:0.5;" />\n'.format(letter_begin+itequery*letter_spacing, subject_spacing+iteline*line_spacing-10))
                                    outf.write('<text x="{}" y="{}" style="font-size:12px;">{}</text>\n'.format(letter_begin+itequery*letter_spacing, subject_spacing+iteline*line_spacing, tali[itequery]))
                            # <breaker in the alignment PGDNS>
                            else:
                                # <SVG::query>
                                outf.write('<rect x="{}" y="{}" width="7" height="10" style="fill:blue; fill-opacity:0.2;" />\n'.format(letter_begin+itequery*letter_spacing, yquery+iteline*line_spacing-10))
                                outf.write('<text x="{}" y="{}" style="font-size:12px;">{}</text>\n'.format(letter_begin+itequery*letter_spacing, yqueryiteline*+line_spacing, aa))
                                # <SVG::subject>
                                outf.write('<text x="{}" y="{}" style="font-size:12px;">{}</text>\n'.format(letter_begin+itequery*letter_spacing, subject_spacing+iteline*line_spacing, tali[itequery]))
                        else:
                            # <SVG::query>
                            outf.write('<text x="{}" y="{}" style="font-size:12px;">{}</text>\n'.format(letter_begin+itequery*letter_spacing, yquery+iteline*line_spacing, aa))
                            # <SVG::subject>
                            outf.write('<text x="{}" y="{}" style="font-size:12px;">{}</text>\n'.format(letter_begin+itequery*letter_spacing, subject_spacing+iteline*line_spacing, tali[itequery]))
                        
                outf.write("</svg>")
    
### I/O
    
def read_fasta(query, proteins, ffdata, ffindex): 
    """ get protein sizes and fasta sequences
    """
    # shorter name for ffindex
    shorter_prots = dict()
    for prot in proteins:
        shorter_prots[prot.split("|")[1]] = prot
    # read ffindex
    data_ffindex = dict()
    with open(ffindex) as inf:
        for line in inf:
            tmp = line.split()
            if tmp[0] in shorter_prots:
                name = shorter_prots[tmp[0]]
                data_ffindex[name] = (int(tmp[1]), int(tmp[2]))

    # get seq and sizes
    sizes = dict()
    fasta = dict()
    with open(ffdata, "rb") as inf:
        for prot in data_ffindex:
            start, size = data_ffindex[prot]
            inf.seek(start)
            seqpart = inf.read(size).decode()
            tokeep = False
            for line in seqpart.split("\n"):
                if line == "\x00":
                    continue
                if line[0] == ">":
                    name = line[1:].split()[0]
                    tokeep = False
                    if name == prot:
                        tokeep = True
                        sizes[name] = 0
                        fasta[name] = ""
                elif tokeep:
                    fasta[name] += line.strip()
                    sizes[name] += len(line.strip())
    # read query
    with open(query) as inf:
        cnt = 0
        for line in inf:
            if line[0] == ">" and cnt == 0:
                fasta["query"] = ""
                sizes["query"] = 0
                cnt += 1
            elif line[0] == ">" and cnt > 0:
                print ("Error, multiple fasta sequences in query {}".format(query), file=sys.stderr)
                sys.exit(1)
            else:
                fasta["query"] += line.strip()
                sizes["query"] += len(line.strip())
    return fasta, sizes

def local_fasta(query, proteins, ffdata, ffindex):
    fasta, sizes = dict(), dict()
    #query sequence
    with open(query) as inf:
        cnt = 0
        for line in inf:
            if line[0] == ">" and cnt == 0:
                fasta["query"] = ""
                sizes["query"] = 0
                cnt += 1
            elif line[0] == ">" and cnt > 0:
                print ("Error, multiple fasta sequences in query {}".format(query), file=sys.stderr)
                sys.exit(1)
            else:
                fasta["query"] += line.strip()
                sizes["query"] += len(line.strip())
    # other sequences
    for prot in proteins:
        size = np.random.randint(100, 300)
        seq = "".join([np.random.choice(AA1) for i in range(size)])
        fasta[prot] = seq
        sizes[prot] = size
    return fasta, sizes

def writeCSS(workdir):
    """
    brief write CSS style
    """
    
    path = os.path.join(workdir,"style.css")
    with open(path, "w") as outf:
        outf.write(""".caption{
            font-size: 20px;
            font-family: Georgia, serif;
        }
        .tooltip{
            font-size: 15px;
        }
        .tooltip_bg{
            fill: white;
            stroke: black;
            stroke-width: 1.2;
            opacity: 0.85;
        } 
        """)

def writeAjax(workdir):
    """
    """
    
    path = os.path.join(workdir,"script.js")
    with open(path, "w") as outf:
        outf .write("""function init(evt) {
            if ( window.svgDocument == null ) {
                svgDocument = evt.target.ownerDocument;
            }
            tooltip = svgDocument.getElementById('tooltip');
            tooltip_bg = svgDocument.getElementById('tooltip_bg');
        }
    
        function ShowTooltip(evt, mouseovertext){
            tooltip.setAttributeNS(null,"x",evt.pageX+11);
            tooltip.setAttributeNS(null,"y",evt.pageY+13);
            tooltip.firstChild.data = mouseovertext;
            tooltip.setAttributeNS(null,"visibility","visible");
    
            length = tooltip.getComputedTextLength();
            tooltip_bg.setAttributeNS(null,"width",length+8);
            tooltip_bg.setAttributeNS(null,"x",evt.pageX+8);
            tooltip_bg.setAttributeNS(null,"y",evt.pageY+0);
            tooltip_bg.setAttributeNS(null,"visibility","visible");
        }
    
        function HideTooltip(evt) {
            tooltip.setAttributeNS(null,"visibility","hidden");
            tooltip_bg.setAttributeNS(null,"visibility","hidden");
        }
        """)

def write_html(order, data, output, workdir):
    """ write final html file
    """
    dirfile = os.path.dirname(output)
    fullwork = os.path.abspath(workdir)
    relpath = os.path.relpath(fullwork, dirfile)
    relpath = os.path.join(relpath, "svg")
    
    with open(output, "w") as outf:
        itehit = 0
        outf.write(html_header())
        outf.write('        <ol>\n')
        prevprot = order[0][0]
        for prot, hitnum in order:
            if hitnum == 0:
                # write domain arrangement of previous prot
                if prevprot != prot:
                    outf.write('            <embed src="{}" name="svgmap" type="image/svg+xml">\n'.format(os.path.join(relpath, prevprot+".svg")))
                prevprot = prot
                # write protein header
                headerprot = create_protheader(prot)
                outf.write(headerprot)
            # start hit
            outf.write('            <ul>\n')
            # hit score
            outf.write("                <li><b>Eval={} ; Probab={} ; Id={} ; Sim={} ; </b></br>\n".format(data[prot][hitnum]["E-value"], data[prot][hitnum]["Probab"], data[prot][hitnum]["Identities"], data[prot][hitnum]["Similarity"]))
            # alignment
            outf.write('                <a href="#{0}" onclick="showhide(\'div{0}\', \'{1}\');">Show/hide Alignment</a>\n'.format(itehit, itehit, os.path.join(relpath, "{}_{}.svg".format(prot, hitnum))))
            outf.write('                <div id="div{0}" style="display: none;" name="{0}"></div>\n'.format(itehit))
            # HCA 
            outf.write('                <a href="#{0}" onclick="showhide(\'hca_s{0}\', \'{1}\');">HCA Subject</a>\n'.format(itehit, os.path.join(relpath, "hca_t_{}_{}.svg".format(prot, hitnum))))
            outf.write('                <div id="hca_s{0}" style="display: none;" name="hca_s{0}"></div>\n'.format(itehit))
            outf.write('                <a href="#{0}" onclick="showhide(\'hca_q{0}\', \'{1}\');">HCA Query</a>\n'.format(itehit, os.path.join(relpath, "hca_q_{}_{}.svg".format(prot, hitnum))))
            outf.write('                <div id="hca_q{0}" style="display: none;" name="hca_q{0}"></div>\n'.format(itehit))
            # end of hit
            outf.write('            </ul>\n')
            itehit += 1
        # last protein
        outf.write('            <embed src="{}" name="svgmap" type="image/svg+xml">\n'.format(os.path.join(relpath, prot+".svg")))    
        outf.write('        </ol>\n')
        # html footer
        outf.write(html_footer())
        
def main():
    # get parameters
    params = get_cmd()
    
    if not os.path.isdir(params.workdir):
        os.makedirs(params.workdir)
    # normalize outfile path and outdir path
    fullwork = os.path.abspath(params.workdir)
    params.workdir = fullwork
            
    # read tremolo result
    order, hits, domains, proteins = read_tremolo(params.inputres, fetch="all")
    
    # get sequences and sizes
    fasta, sizes = read_fasta(params.query, proteins, params.ffdata, params.ffindex)
    #fasta, sizes = local_fasta(params.query, proteins, params.ffdata, params.ffindex)
    
    # compute hca of sequences
    dhca = dict()
    for prot in fasta:
        dhca[prot] = segmentation(fasta[prot])
    
    # create domain arrangements and alignments svg
    svg_tremolo(fasta, hits, domains, dhca, sizes, params.workdir)
                        
    # finalize html
    write_html(order, hits, params.output, params.workdir)

    sys.exit(0)

if __name__ == "__main__":
    main()

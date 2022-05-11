import os, sys
import random
import csv
from numpy import random, argmax, diff, mean
import plotly
import plotly.graph_objs as go
from plotly import tools
# from plotly.subplots import make_subplots

# reload(sys)  
# sys.setdefaultencoding('utf8')

# Check the command line parameters
if len(sys.argv) < 4:
    print ("One or more parameters are missing!")
    sys.exit(1)

global prefix
global chromosmes
global FileFolder
global out_folder

prefix, chromo, FileFolder, out_folder = sys.argv[1:5]    

chromosmes = chromo.split(',')

def ListChanger(t):
    StartingBool = [True]*len(chromosmes) + [False]*(NEvents) + [True]
    for i in t:
        StartingBool[i] = True  
    return StartingBool  

fig=tools.make_subplots(rows=1, cols=1, shared_xaxes=False,
                          shared_yaxes=False, vertical_spacing=0.1, print_grid=False)
# fig= make_subplots(rows=1, cols=1, shared_xaxes=False,
#                           shared_yaxes=False, vertical_spacing=0.1, print_grid=False)



cols=['rgba(166,206,227,0.4)','rgba(31,120,180,0.4)','rgba(178,223,138,0.4)','rgba(51,160,44,0.4)','rgba(251,154,153,0.4)','rgba(227,26,28,0.4)','rgba(253,191,111,0.4)','rgba(255,127,0,0.4)','rgba(202,178,214,0.4)','rgba(106,61,154,0.4)','rgba(204,204,0,0.4)']
last_pos=0
last_idx=len(chromosmes)
c_counter=0
difference=0
mean_val_for_axes=0
plotAnnotesy1 = []
FilteringIndeces = [[] for i in range(10)]
seg_plot = []
supp_data = []
global NEvents
NEvents = 0
for idx,c in enumerate(chromosmes):
    Flag=False
    DataTab = []
    Vcf = os.path.join(FileFolder, str(c) + '_HSLMResults_' + prefix + '.txt')
    CallFile = os.path.join(FileFolder, str(c) + '_FractionCallResults_' + prefix + '.txt')
    with open(Vcf, 'r') as f:
        # reader = csv.reader(f, delimiter='\t')
        reader = csv.reader(f, delimiter='\t')
        next(f)
        for row in reader:
            DataTab.append([ row[i] for i in [0, 1, 4, 5] ])
    CallsTab = []
    with open(CallFile, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row[0] == ('V1'):    
                CallsTab.append([ row[i] for i in [1, 2, 4, 5, 6, 7, 8] ])
    NEvents += len(CallsTab)
    # pos_vec=map(int, zip(*DataTab)[1])
    # pos_vec=map(int, list(zip(*DataTab))[1])
    pos_vec=list(map(int, list(zip(*DataTab))[1]))
    if idx > 0:
        difference=(pos_vec[0] + last_pos)-(last_pos+1)
    # seg=map(float,zip(*DataTab)[2])
    
    # seg=map(float,list(zip(*DataTab))[2])
    seg=list(map(float,list(zip(*DataTab))[2]))

    if mean([mean_val_for_axes, max(seg)]) > mean_val_for_axes:
        # mean_val_for_axes=mean(seg)
        mean_val_for_axes=mean(seg)
        mean_val_for_axes=mean([mean_val_for_axes, max(seg)])
    new_pos_vec=[i+last_pos-difference for i in pos_vec]
    # seg_plot = seg_plot + zip(new_pos_vec, zip(*DataTab)[3])
    seg_plot = seg_plot + list(zip(new_pos_vec, list(zip(*DataTab))[3]))
    last_pos=new_pos_vec[-1]+7000000
    GapIdx=argmax(diff(new_pos_vec))
    rep_idxs=[GapIdx, GapIdx+1]
    # new_pos_dict = dict(zip(pos_vec, new_pos_vec))
    new_pos_dict = dict(list(zip(pos_vec, new_pos_vec)))
    repl=[None, None]
    for repidx, rl in list(zip(rep_idxs, repl)):
        new_pos_vec[repidx]=rl
        seg[repidx]=rl
    if idx in [10, 20, 30, 40]:
       c_counter=0
    trace5=go.Scatter(
        x=new_pos_vec,
        y=[round(i,3) for i in seg if i is not None],
        name=str(c),
        mode='lines',
        showlegend=False,
        text=list(map(lambda x:str(x), [round(i,2) for i in seg if i is not None])),
        hoverinfo="x+text+name",
        line=dict(width=0.5,
            color=cols[c_counter])       
    )
    
    fig.append_trace(trace5, 1, 1)
    if len(DataTab[0][0])>2:
        ChrOut = DataTab[0][0][3:]
    else:
        ChrOut = DataTab[0][0]
    
    for ind, ele in enumerate(CallsTab):
        if ele[2] in ['Deletion', 'Double Deletion']:fcol = 'rgba(139, 0, 0, 0.4)'; textpos = 'bottom'
        else: fcol = 'rgba(50, 171, 96, 0.4)'; textpos = 'top'
        calls2BAF = go.Scatter(x=[new_pos_dict.get(int(ele[0])), new_pos_dict.get(int(ele[1]))],y=[ele[6], ele[6]], mode='lines', hoverinfo='skip', fill = 'tozeroy', showlegend=False, line=dict(width=1, color=fcol))
        Idx = int(round(float(ele[4]),1)*10)
        
        if Idx > 0:
            Idx = Idx-1
        if Idx >= 10:
            Idx = 9

        FilteringIndeces[Idx].append(ind+last_idx)
        supp_data.append(calls2BAF)
        Linkline = '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr' + ChrOut + ':' + str(ele[0]) + '-' + str(ele[1]) + '">{}</a>' 
        plotAnnotesy1.append(dict(x = (int(new_pos_dict.get(int(ele[0]))) + int(new_pos_dict.get(int(ele[1]))))/2,
                             y = float(ele[6])+0.05,
                             text = Linkline.format('<b>' + u'\u24d8' + '</b>'),
                             showarrow = True,
                             xanchor = 'center',
                             yanchor = 'top',
                             visible = True,
                             font = dict(color = 'rgb(0,0,0,0)'),
                             hovertext = "<b>POS:</b> %s:%s-%s <br> <b>CN:</b> %s <br> <b>AF:</b> %s &plusmn; %s" % (str(ChrOut), ele[0], ele[1], ele[3], ele[4], ele[5]),
                             hoverlabel = dict(bgcolor=fcol)
                             ))
    c_counter += 1
    last_idx = len(list(zip(*CallsTab))[0])+last_idx
for sd in supp_data:
    fig.data.append(sd)
plotAnnotesy1.append(dict(text='Allele Frequency', x=0, y=-0.11, showarrow=False, align='left', yref='paper', xref = 'paper', xanchor = 'left', yanchor = 'bottom', visible = True))

for l in range(len(FilteringIndeces)):
    FilteringIndeces[l] = list(set(FilteringIndeces[l]))

trace6=go.Scatter(
    x=list(map(int, list(zip(*seg_plot))[0])),
    y=[round(i,2) for i in list(map(float, list(zip(*seg_plot))[1])) if i is not None],
    name='Segmentation',
    mode='lines',
    showlegend=True,
    text=list(map(lambda x:str(x), [round(i,3) for i in list(map(float, list(zip(*seg_plot))[1])) if i is not None])),
    hoverinfo="text+name",
    line=dict(width=1, color='red')   
)
fig.data.append(trace6)

updatemenus=list([
        dict(
            buttons=list([
                dict(
                    args=[{'visible': ListChanger([item for sublist in FilteringIndeces for item in sublist])}, {'annotations': list(map(lambda x:plotAnnotesy1[x], [item-len(chromosmes) for sublist in FilteringIndeces for item in sublist] + [len(plotAnnotesy1)-1]))}],
                    label='All',
                    method='update'
                ),
                dict(
                    args=[{'visible': ListChanger([item for sublist in FilteringIndeces[1:] for item in sublist])}, {'annotations': list(map(lambda x:plotAnnotesy1[x], [item-len(chromosmes) for sublist in FilteringIndeces[1:] for item in sublist] + 
[len(plotAnnotesy1)-1]))}],
                    label='10%',
                    method='update'
                ),
                dict(
                    args=[{'visible': ListChanger([item for sublist in FilteringIndeces[2:] for item in sublist])}, {'annotations': list(map(lambda x:plotAnnotesy1[x], [item-len(chromosmes) for sublist in FilteringIndeces[2:] for item in sublist] + 
[len(plotAnnotesy1)-1]))}],
                    label='20%',
                    method='update'
                ),
                dict(
                    args=[{'visible': ListChanger([item for sublist in FilteringIndeces[3:] for item in sublist])}, {'annotations': list(map(lambda x:plotAnnotesy1[x], [item-len(chromosmes) for sublist in FilteringIndeces[3:] for item in sublist] + 
[len(plotAnnotesy1)-1]))}],
                    label='30%',
                    method='update'
                ),
                dict(
                    args=[{'visible': ListChanger([item for sublist in FilteringIndeces[4:] for item in sublist])}, {'annotations': list(map(lambda x:plotAnnotesy1[x], [item-len(chromosmes) for sublist in FilteringIndeces[4:] for item in sublist] + 
[len(plotAnnotesy1)-1]))}],
                    label='40%',
                    method='update'
                ),
                dict(
                    args=[{'visible': ListChanger([item for sublist in FilteringIndeces[5:] for item in sublist])}, {'annotations': list(map(lambda x:plotAnnotesy1[x], [item-len(chromosmes) for sublist in FilteringIndeces[5:] for item in sublist] + 
[len(plotAnnotesy1)-1]))}],
                    label='50%',
                    method='update'
                ),
                dict(
                    args=[{'visible': ListChanger([item for sublist in FilteringIndeces[6:] for item in sublist])}, {'annotations': list(map(lambda x:plotAnnotesy1[x], [item-len(chromosmes) for sublist in FilteringIndeces[6:] for item in sublist] + 
[len(plotAnnotesy1)-1]))}],
                    label='60%',
                    method='update'
                ),  
                dict(
                    args=[{'visible': ListChanger([item for sublist in FilteringIndeces[7:] for item in sublist])}, {'annotations': list(map(lambda x:plotAnnotesy1[x], [item-len(chromosmes) for sublist in FilteringIndeces[7:] for item in sublist] + 
[len(plotAnnotesy1)-1]))}],
                    label='70%',
                    method='update'
                ),
                dict(
                    args=[{'visible': ListChanger([item for sublist in FilteringIndeces[8:] for item in sublist])}, {'annotations': list(map(lambda x:plotAnnotesy1[x], [item-len(chromosmes) for sublist in FilteringIndeces[8:] for item in sublist] + 
[len(plotAnnotesy1)-1]))}],
                    label='80%',
                    method='update'
                ),
                dict(
                    args=[{'visible': ListChanger([item for sublist in FilteringIndeces[9:] for item in sublist])}, {'annotations': list(map(lambda x:plotAnnotesy1[x], [item-len(chromosmes) for sublist in FilteringIndeces[9:] for item in sublist] + 
[len(plotAnnotesy1)-1]))}],
                    label='90%',
                    method='update'
                ), 
                dict(
                    args=[{'visible': ListChanger([item for sublist in [] for item in sublist])}, {'annotations': list(map(lambda x:plotAnnotesy1[x], [item-len(chromosmes) for sublist in [] for item in sublist] + [len(plotAnnotesy1)-1]))}],
                    label='None',
                    method='update'
                ),                      
            ]),
            direction = 'left',
            font = dict(size=10),
            pad = {'l':0, 'b': 0},
            showactive = True,
            type = 'buttons',
            x = 0,
            y = -0.16,
            xanchor = 'left',
            yanchor = 'bottom'
        ),
    ])


fig['layout'].update(
    title='Calling Results',
    updatemenus=updatemenus,
    annotations = plotAnnotesy1,
    xaxis=dict(showspikes=False, title='Position (bp)', range=[0, last_pos+500000], exponentformat='SI'),
    yaxis=dict(showspikes=False, title='log2ratio', zeroline=True),
    margin=dict(l=100,r=100,t=100,b=100),
    hovermode='closest',
)


config={'displaylogo': False,
          'modeBarButtonsToRemove': ['toImage', 'sendDataToCloud', 'select2d', 'lasso2d',
                                     'hoverClosestCartesian', 'hoverCompareCartesian','toggleHover', 'zoomIn2d', 'zoomOut2d', 'toggleSpikelines']}

out_stats = os.path.join(out_folder, (prefix + '_CNV_calling_results.html'))
plotly.offline.plot(fig, filename=out_stats, auto_open=False, config=config, show_link=False)

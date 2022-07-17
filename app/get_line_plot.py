import os
import glob
import numpy as np
import pandas as pd
import json
import plotly.graph_objs as go
import plotly.offline as plotly
import plotly.figure_factory as ff
import math

def plot(xx, yy, name,color):
    print("x",xx)
    print("y",yy)
    plot = go.Scatter(
        x=xx,
        y=yy,
        
        mode = "lines",
        # marker=dict(color=color,size=30),
        line={'color': color, 'width': 5},
        # line=dict(color=color),
        name=name
        # marker=dict(color=colors[key])
    )

    return plot

############################################
filenbr="8"
check_content=1 # 0: angle; 1: kn; 2:kg
color_list=['#00b894','#0984e3','#d63031']
target=""
yname=""
yrange=[]
ynticks=5
if check_content==0:
    target="angle"
    yrange=[0,90]
    yname=target
    ynticks=5
if check_content==1:
    target="kn"
    yrange=[0.18,0.22]
    yname=target
    ynticks=10
if check_content==2:
    target="kg"
    yrange=[0,0.5]
    yname=target
    ynticks=10

##########################################33


filepath ='../build/'


filename=''
counter=0
plots=[]
x1=[]
y1=[]
y2=[]
y3=[]
filename=filepath+filenbr+"_info.csv"
tdf=pd.read_csv(filename)
x1=[*range(1,21)]
y1=tdf[target].values

y1=np.absolute(y1);
plots.append(plot(x1, y1, target,color_list[0]))


layout= go.Layout(
    legend=dict(x=0.5, y=1),
    font=dict(
        size=36
    ),
    xaxis=dict(
        title="the points on the line",
        range=[0,20],
        # exponentformat='power',
        showticksuffix='all',
        showtickprefix='all',
        showexponent='all',
        # autotick=True,
        nticks=20,
        tickfont=dict(
            size=20
        ),
        # type='log',
        showline=True, 
        linewidth=2, 
        showgrid=True, 
        gridwidth=1, 
        gridcolor='grey',

        linecolor='black'
        # ticktext=xticks
    ),
    yaxis=dict(
        title=yname,
        # tickformat='.1e',
        # exponentformat='power',
        # ticks='',
        # range=[0,5000],
        # tick0=0,
        # dtick=1,
        # tickangle=-45,
        tickfont=dict(
            size=20
        ),
        # autorange=True,
        showline=True, 
        linewidth=2, 
        linecolor='black',
        showgrid=True, 
        gridwidth=1, 
        gridcolor='grey',
        # tick0=1, 
        # dtick=1,
        # nticks=10,
        range=yrange,
        #############################################
        
        
        nticks=ynticks,
        
        ###############################################
    ),
    # font=dict(
    #     size=16
    # ),
    hovermode='closest',
    paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)'
)




fig = go.Figure(data=plots,layout=layout)
# if check_content==2 or check_content==3:
#     fig.update_yaxes(
#     # type='log'
#     )
# fig.update_layout(yaxis_tickformat = '%')
fig.write_image(filenbr+"_info_"+target + ".svg")


# fig.show()
#if output is not None:
#    fig.write_image(output + ".svg")
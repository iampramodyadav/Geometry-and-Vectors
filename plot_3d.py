"""
@Author:    Pramod Kumar Yadav
@email:     pkyadav01234@gmail.com
@Date:      August, 2023
@status:    development
@PythonVersion: python3

"""
# Import dependencies
import plotly
import plotly.graph_objs as go
import numpy as np
import pandas as pd
#---------------------------------------

def plot_3d_point(list):
    """
    Plots a 3D point cloud.

    Args:
        list (list): A list of 3D points, each point is a list of three numbers.

    Returns:
        None.
    """    
    x_data=[list[i][0] for i in range(len(list))]
    y_data=[list[i][1] for i in range(len(list))]
    z_data=[list[i][2] for i in range(len(list))]
    
    trace = go.Scatter3d(
        x=x_data,
        y=y_data,
        z=z_data,
        mode ='markers',
        marker = dict(size= 10,opacity= 0.9,color=z_data, colorscale='plotly3'))
    
    # Configure the layout.
    layout = go.Layout(margin={'l': 0, 'r': 0, 'b': 0, 't': 10})
    data = [trace]
    plot_figure = go.Figure(data=data, layout=layout)
    plot_figure.update_layout(width=600,height=500,)
    # plot_figure.show()
    return plot_figure

#---------------------------------------
def plot_arrow_tip(point_pair):
    """
    This function plots a vector tip cone from head and tail.

    Args:
        point_pair (list): A list of two points, the first point is the tail and the second point is the head.

    Returns:
        fig (go.Figure): A plotly figure object.
    """
    x_data=[point_pair[i][0] for i in range(len(point_pair))]
    y_data=[point_pair[i][1] for i in range(len(point_pair))]
    z_data=[point_pair[i][2] for i in range(len(point_pair))]

    vec_cone = np.array(point_pair[1])-np.array(point_pair[0])
    # print(vec_cone)
    vec_norm = np.linalg.norm(vec_cone)
    
    u = vec_cone[0]/vec_norm
    v = vec_cone[1]/vec_norm
    w = vec_cone[2]/vec_norm

    fig = go.Figure(data=go.Cone(
        x=[x_data[1]],
        y=[y_data[1]],
        z=[z_data[1]],
        
        u=[u],
        v=[v],
        w=[w],
        
        showscale=False,
        sizemode="scaled",
        sizeref=0.1,
        anchor="tip"))
    
    fig.update_layout(title = '3D Vector Tip Plot',width=500,height=400,)
    
    # fig.show()
    return fig
    
def plot_3d_line(list):
    """
    Plots a 3D line.

    Args:
        list (list): A list of 3D points, each point is a list of three numbers.

    Returns:
        go.Figure.
    """    
    x_data=[list[i][0] for i in range(len(list))]
    y_data=[list[i][1] for i in range(len(list))]
    z_data=[list[i][2] for i in range(len(list))]
    
    trace=go.Scatter3d(
        x=x_data, 
        y=y_data, 
        z=z_data,
        # mode='markers',
        marker=dict(size=5,
                    opacity= 0.5,
                    # color=z_data,
                    # colorscale='Viridis',
                   ),
        line=dict(color='darkblue',width=2)
    )
    # Configure the layout.
    layout = go.Layout(margin={'l': 0, 'r': 0, 'b': 0, 't': 30})
    data = [trace]
    plot_figure = go.Figure(data=data, layout=layout)
    plot_figure.update_layout(title = '3D Line Plot',width=500,height=400,)
    # plot_figure.show()
    return plot_figure
    
#---------------------------------------

def plot_lines_from_points(first_pair, *list_pair):
    """
    Plots a series of 3D lines from a list of 3D points.

    Args:
        list_pair (list): lists of 3D points. Each list is a list of points that form a line.

    Returns:
        go.Figure.
    """
    fig1  = plot_3d_line(first_pair)
    fig0 = plot_arrow_tip(first_pair)
    fig = go.Figure(data=fig1.data + fig0.data)
    
    for pair in list_pair:
        fig1= plot_3d_line(pair)
        fig0 = plot_arrow_tip(pair)
        
        fig2 = go.Figure(data=fig.data + fig1.data + fig0.data)
        fig = fig2

    fig.update_layout(margin={'l': 0, 'r': 0, 'b': 0, 't': 30},title = '3D Vector Plot',width=700,height=400,)
    return fig

#---------------------------------------

def surf_plot(x_data, y_data, z_data):
    """
    Plots a 3D surface plot.

    Args:
        x_data (list): A list of x-coordinates.
        y_data (list): A list of y-coordinates.
        z_data (list): A list of z-coordinates.
    Returns:
        go.Figure.
    """
    
    trace = go.Surface(x = x_data, y = y_data, z =z_data )    
    # Configure the layout.
    data = [trace]
    layout = go.Layout(title = '3D Surface plot',width=700,height=400,margin={'l': 0, 'r': 0, 'b': 0, 't': 30})
    fig = go.Figure(data = data, layout=layout)
    fig.show()

if __name__ == "__main__":
    #------------------------------
    list_point=[[1,3,3],
          [4,3,6],
          [7,5,9],
          [10,11,13]]
    fig1=plot_3d_point(list_point)
    fig1.show()
    #------------------------------
    list1=[[0,0,0],[2,1,1]]
    fig0=plot_arrow_tip(list1)
    fig0.show()
    #------------------------------
    list=[[1,2,3],[4,6,6],[7,8,9],[10,11,12],[2,5,6]]
    fig2=plot_3d_line(list)
    fig2.show()
    #------------------------------
    list1=[[0,0,0],[1,0,0]]
    list2=[[0,0,0],[0,1,0]]   
    list3=[[0,0,0],[0,0,1]] 
    list4=[[0,0,0],[1,1,1]]
    fig3=plot_lines_from_points(list1,list2,list3,list4)
    fig3.show()





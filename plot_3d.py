# Import dependencies
import plotly
import plotly.graph_objs as go

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
    plot_figure.update_layout(width=800,height=700,)
    plot_figure.show()

#---------------------------------------

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
        marker=dict(size=4,opacity= 0.9,color=z_data,colorscale='Viridis',),
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
        list_pair (list): A list of lists of 3D points. Each sublist is a list of points that form a line.

    Returns:
        go.Figure.
    """
    fig= plot_3d_line(first_pair)
    
    for pair in list_pair:
        fig1= plot_3d_line(pair)
        fig2 = go.Figure(data=fig.data + fig1.data)
        fig = fig2

    fig.update_layout(margin={'l': 0, 'r': 0, 'b': 0, 't': 30},title = '3D Lines Plot',width=800,height=400,)
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
    layout = go.Layout(title = '3D Surface plot',width=500,height=400,margin={'l': 0, 'r': 0, 'b': 0, 't': 30})
    fig = go.Figure(data = data, layout=layout)
    fig.show()







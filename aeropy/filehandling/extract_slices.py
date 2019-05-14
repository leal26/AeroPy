import numpy as np
# import the simple module from the paraview
from paraview.simple import *
# disable automatic camera reset on 'Show'
# paraview.simple._DisableFirstRenderCameraReset()

# create a new 'STL Reader'
directory = 'D:\\GitHub\\AeroPy\\examples\\JAXA_files'
filename = 'wing_right'
input = directory + '\\raw\\' + filename + '.stl'
output = directory + '\\processed\\' + filename + '_'
fuselagestl = STLReader(FileNames=[input])

# get active view
# renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1256, 816]

# show data in view
# fuselagestlDisplay = Show(fuselagestl, renderView1)

# reset view to fit data
# renderView1.ResetCamera()

# show color bar/color legend
# fuselagestlDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'STLSolidLabeling'
# sTLSolidLabelingLUT = GetColorTransferFunction('STLSolidLabeling')

# get opacity transfer function/opacity map for 'STLSolidLabeling'
# sTLSolidLabelingPWF = GetOpacityTransferFunction('STLSolidLabeling')

# create a new 'Transform'
transform1 = Transform(Input=fuselagestl)

# Properties modified on transform1.Transform
transform1.Transform.Rotate = [0.0, -2.3067, 0.0]


def cosine_spacing(start, stop, num=50, offset=0):
    # calculates the cosine spacing
    index = np.linspace(0., 1., num)
    spacing = .5*(1.-np.cos(np.pi*(index-offset)))

    points = start+spacing*(stop-start)

    return points


# x_slice = cosine_spacing(0.0001, 38.65, 100)
x_slice = cosine_spacing(0.0001, 38.6736, 100)
x_slice = np.linspace(0.8, 4.57, 100)

for i, x in enumerate(x_slice):
    # create a new 'Slice'
    slice1 = Slice(Input=transform1)

    # Properties modified on slice1.SliceType
    slice1.SliceType.Origin = [x, 0.0, 0.0]

    # show data in view
    # slice1Display = Show(slice1, renderView1)

    # hide data in view
    # Hide(fuselagestl, renderView1)

    # show color bar/color legend
    # slice1Display.SetScalarBarVisibility(renderView1, True)

    # save data
    SaveData(output+str(i)+'.csv', proxy=slice1, Precision=12)

# saving camera placements for all active views

# current camera placement for renderView1
# renderView1.CameraPosition = [19.325271606445312, -1.9311904907226562e-05, 73.83615578949325]
# renderView1.CameraFocalPoint = [19.325271606445312, -1.9311904907226562e-05, -1.0328417336568236]
# renderView1.CameraParallelScale = 19.37752244672469

# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).

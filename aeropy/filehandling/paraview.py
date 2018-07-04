#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

network_number = 2
filename = 'test_network'
directory = 'C:\\Users\\leal26\\Documents\\GitHub\\AeroPy\\aeropy\\CST\\'

# get active view
renderView = GetActiveViewOrCreate('RenderView')

assembly = []
for i in range(1,network_number+1):
    # create a new 'XML Structured Grid Reader'
    test_network_vts = XMLStructuredGridReader(FileName=[directory + filename + str(i)+'.vts'])

    # show data in view
    test_network_vtsDisplay = Show(test_network_vts, renderView)
    # trace defaults for the display properties.
    test_network_vtsDisplay.Representation = 'Surface With Edges'
    test_network_vtsDisplay.ColorArrayName = [None, '']
    test_network_vtsDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    test_network_vtsDisplay.SelectOrientationVectors = 'None'
    test_network_vtsDisplay.ScaleFactor = 0.1
    test_network_vtsDisplay.SelectScaleArray = 'None'
    test_network_vtsDisplay.GlyphType = 'Arrow'
    test_network_vtsDisplay.GlyphTableIndexArray = 'None'
    test_network_vtsDisplay.DataAxesGrid = 'GridAxesRepresentation'
    test_network_vtsDisplay.PolarAxes = 'PolarAxesRepresentation'
    test_network_vtsDisplay.ScalarOpacityUnitDistance = 0.3272506722223079

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    test_network_vtsDisplay.OSPRayScaleFunction.Points = [2.326428429822192, 0.0, 0.5, 0.0, 37.626781425423815, 1.0, 0.5, 0.0]

# reset view to fit data
renderView.ResetCamera()

# update the view to ensure updated data information
renderView.Update()

#### saving camera placements for all active views

# current camera placement for renderView1
renderView.CameraPosition = [0.12476075744808501, 3.1845058646858693, 0.3710215545807592]
renderView.CameraFocalPoint = [0.5, 0.5, 0.0037752263491506906]
renderView.CameraViewUp = [-0.30729811760225784, -0.17101732138568032, 0.9361201539888863]
renderView.CameraParallelScale = 0.7079657120931511

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
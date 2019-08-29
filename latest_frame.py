
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XDMF Reader'
outputxdmf = XDMFReader(FileNames=['./output.xdmf'])
outputxdmf.PointArrayStatus = ['f_5-0']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# set active source
SetActiveSource(outputxdmf)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [859, 567]

# show data in view
outputxdmfDisplay = Show(outputxdmf, renderView1)

# get color transfer function/color map for 'f_50'
f_50LUT = GetColorTransferFunction('f_50')
f_50LUT.RGBPoints = [0.48044341802597046, 0.231373, 0.298039, 0.752941, 0.4993125796318054, 0.865003, 0.865003, 0.865003, 0.5181817412376404, 0.705882, 0.0156863, 0.14902]
f_50LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'f_50'
f_50PWF = GetOpacityTransferFunction('f_50')
f_50PWF.Points = [0.48044341802597046, 0.0, 0.5, 0.0, 0.5181817412376404, 1.0, 0.5, 0.0]
f_50PWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
outputxdmfDisplay.Representation = 'Surface'
outputxdmfDisplay.ColorArrayName = ['POINTS', 'f_5-0']
outputxdmfDisplay.LookupTable = f_50LUT
outputxdmfDisplay.OSPRayScaleArray = 'f_5-0'
outputxdmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
outputxdmfDisplay.SelectOrientationVectors = 'None'
outputxdmfDisplay.ScaleFactor = 4.0
outputxdmfDisplay.SelectScaleArray = 'f_5-0'
outputxdmfDisplay.GlyphType = 'Arrow'
outputxdmfDisplay.GlyphTableIndexArray = 'f_5-0'
outputxdmfDisplay.GaussianRadius = 0.2
outputxdmfDisplay.SetScaleArray = ['POINTS', 'f_5-0']
outputxdmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
outputxdmfDisplay.OpacityArray = ['POINTS', 'f_5-0']
outputxdmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
outputxdmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
outputxdmfDisplay.SelectionCellLabelFontFile = ''
outputxdmfDisplay.SelectionPointLabelFontFile = ''
outputxdmfDisplay.PolarAxes = 'PolarAxesRepresentation'
outputxdmfDisplay.ScalarOpacityFunction = f_50PWF
outputxdmfDisplay.ScalarOpacityUnitDistance = 2.418271175121957

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
outputxdmfDisplay.DataAxesGrid.XTitleFontFile = ''
outputxdmfDisplay.DataAxesGrid.YTitleFontFile = ''
outputxdmfDisplay.DataAxesGrid.ZTitleFontFile = ''
outputxdmfDisplay.DataAxesGrid.XLabelFontFile = ''
outputxdmfDisplay.DataAxesGrid.YLabelFontFile = ''
outputxdmfDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
outputxdmfDisplay.PolarAxes.PolarAxisTitleFontFile = ''
outputxdmfDisplay.PolarAxes.PolarAxisLabelFontFile = ''
outputxdmfDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
outputxdmfDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show color bar/color legend
outputxdmfDisplay.SetScalarBarVisibility(renderView1, True)

# reset view to fit data
renderView1.ResetCamera()

animationScene1.GoToLast()

# Rescale transfer function
f_50LUT.RescaleTransferFunction(0.0, 1.0)

# Rescale transfer function
f_50PWF.RescaleTransferFunction(0.0, 1.0)

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# get the material library
materialLibrary1 = GetMaterialLibrary()

# hide color bar/color legend
outputxdmfDisplay.SetScalarBarVisibility(renderView1, False)

# current camera placement for renderView1
renderView1.CameraPosition = [20.0, 20.0, 109.2820323027551]
renderView1.CameraFocalPoint = [20.0, 20.0, 0.0]
renderView1.CameraParallelScale = 28.284271247461902

# save screenshot
SaveScreenshot('C:/Users/CE-KPI15/Projects/CHSAFT/latest_frame.png', renderView1, ImageResolution=[1627, 822], 
    # PNG options
    CompressionLevel='0')


import qupath.lib.regions.*
import qupath.imagej.tools.IJTools
import qupath.imagej.gui.IJExtension
import ij.*
IJExtension.getImageJInstance()
def server = getCurrentServer()
def roi = getSelectedROI()
double downsample = 1
def request = RegionRequest.createInstance(server.getPath(), downsample, roi)
def pathImage = IJTools.convertToImagePlus(server, request)
def imp = pathImage.getImage()


// Convert QuPath ROI to ImageJ Roi & add to open image
def roiIJ = IJTools.convertToIJRoi(roi, pathImage)
imp.setRoi(roiIJ)
// Set the line width for the ROI
roiIJ.setStrokeWidth(15)

// Specify the path to the macro file
String macroFilePath = "/Users/polisev/Desktop/auto_scripts/rectangle_3color.ijm"

// Run the macro from the file path
IJ.runMacroFile(macroFilePath);

from cubevis.plots.cubevis_fitsfigure import CubeVisFITSFigure
from cubevis.idealpy.radio import getaxes
from cubevis.idealpy.fits import xyad
from cubevis.idealpy.radio import momentcube, cube_extract
from cubevis.utils import CubeVisArgumentError
import matplotlib.pyplot as mpl
#from matplotlib.mlab import inside_poly
from matplotlib import path
import numpy

class CubeVis(CubeVisFITSFigure):
    """
    The CubeVis Class is a top level class 
    for holding FITS images and cubes and provides
    an interactive platform for data visualization
    """
    def __init__(self, data, hdu=0, figure=None,
                 subplot=(1,1,1), downsample=False,
                 north=False, convention=None,
                 dimensions=[0, 1], slices=[],
                 gauss_width=2, extra_hdus=[],
                 extra_hdu_titles=[],
                 blue_windows=[], red_windows=[],
                 rms_window=[],
                 sigma_crit=3.0, source_distance=140.0,
                 avgfigure=None, pvangle=None,
                 **kwargs):
        """
        For now blue_windows and red_windows are lists of
        two values for blue and redshifted windows
        """
        CubeVisFITSFigure.__init__(self, data, hdu, figure=figure,
                                   subplot=subplot, downsample=downsample,
                                   north=north, convention=convention,
                                   dimensions=dimensions, slices=slices,
                                   **kwargs)
        self.extra_hdus = extra_hdus
        self.extra_hdu_titles = extra_hdu_titles
        self.blue_windows = blue_windows
        self.red_windows = red_windows
        self.rms_window = rms_window
        self.sigma_crit = sigma_crit
        self.avgfigure = avgfigure
        self.pvangle = pvangle
        self.polygon_vert_blue = []
        self.polygon_vert_red = []
        self.blue_polygon_active = False
        self.red_polygon_active = False
        self.source_distance = source_distance
        self.specfigure = mpl.figure(10)
        self.keystrokes = {#key: ('help-text', 'help-prompt')
            'b': ('Start Box-car average',
                  'Click on other edge of box and press e'),
            'e': ('End Box-car average',
                  'First click b to start box'),
            'B': ('Blue-wing Sigma-weighted plot',
                  'Click on the other edge of box and press E'),
            'R': ('Red-wing Sigma-weighted plot',
                  'Click on the other edge of box and press E'),
            'E': ('End blue or red sigma-weighted plot',
                  'First click B or R to start box'),
            't': ('Popup Spectrum',
                  'Pops up spectrum at this position'),
            'p': ('Position-Velocity Start',
                  'PV Cut: Click on other end of line and press v'),
            'v': ('Position-Velocity End',
                  'First start p to start position-velocity cut'),
            'P': ('Position-Velocity by Angle Start',
                  'PV Cut: Click A to enter angle for PV Cut'),
            'A': ('Position-Velocity by Angle End',
                  'First click P to start position-velocity cut by angle'),
            'd': ('Distance measurement start',
                  'Click D to mark end of line for linear distance'),
            'D': ('Distance measurement end',
                  'First click d to start distance measurement'),
            'G': ('Start definition of blue-shifted polygon',
                  'Keep pressing G to define polygon vertices and press c to close'),
            'g': ('Start definition of red-shifted polygon',
                  'Keep pressing G to define polygon vertices and press c to close'),
            'c': ('Close blue or red-shifted polygon',
                  'First press G or g to start blue or red-shifted polygons'),
            'h': ('Print out help prompt for keystrokes',
                  'Print out help prompt for keystrokes')
            }
        self.start_events()
        
        
    def print_keystroke_help(self):
        for key, (help_text, help_prompt) in self.keystrokes.items():
            print("%s: %s" % (key, help_text))

    def on_keypress(self, event):
        if not event.inaxes: return
        if event.key not in self.keystrokes.keys():
            print("Key %s not a recognized keystroke" % event.key)
            self.print_keystroke_help()
        print(event.xdata, event.ydata)
        if event.key == 'b':
            #starts a box for box-averaged plot
            self.box = []
            self.box.append((event.xdata, event.ydata))
            print(self.keystrokes[event.key][1])
            self.avgfigure = None
        if event.key == 'e':
            if not self.box:
                print(self.keystrokes[event.key][1])
                return
            else:
                self.box.append((event.xdata, event.ydata))
                self.draw_box_and_average()
        if event.key == 'B':
            #starts a box for box-averaged blue sigma-weighted plot
            self.box = []
            self.box.append((event.xdata, event.ydata))
            print(self.keystrokes[event.key][1])
            self.blue = True
            self.red = False
        if event.key == 'R':
            #starts a box for box-averaged blue sigma-weighted plot
            self.box = []
            self.box.append((event.xdata, event.ydata))
            print(self.keystrokes[event.key][1])
            self.red = True
            self.blue = False
        if event.key == 'E':
            if not self.box:
                print(self.keystrokes[event.key][1])
                return
            else:
                self.box.append((event.xdata, event.ydata))
                self.draw_box_and_sigma_average()
        if event.key == 't':
            #popup a spectrum
            self.draw_spectrum(event.xdata, event.ydata)
        if event.key == 'p':
            #position-velocity cut
            self.pvcut = []
            self.pvcut.append((event.xdata, event.ydata))
            print(self.keystrokes[event.key][1])
            self.pvfigure = None
        if event.key == 'v':
            if not self.pvcut:
                print(self.keystrokes[event.key][1])
                return
            else:
                self.pvcut.append((event.xdata, event.ydata))
                self.make_posvel_cut()
        if event.key == 'P':
            #print ra, dec for location
            print(event.xdata, event.ydata)
            r = self._wcs.wcs_pix2world([[event.xdata,event.ydata]], 1)
            ra = r[0][0]
            dec = r[0][1]
            print("RA = %s, Dec= %s" % (ra, dec))
            print("RA = %s, Dec = %s" % (sixty(ra*24./360.), sixty(dec)))
        if event.key == 'A':
            print(event.xdata, event.ydata)
            if self.pvangle is None:
                angle = get_input_with_default("Enter angle (reckoned from +ve X-axis) in degrees for PV cut",
                                               90., typeconvert=float)
            else:
                angle = self.pvangle
            self.make_posvel_cut_angle(event.xdata, event.ydata,
                                       angle)
        if event.key == 'd':
            #starts distance cut
            self.distance = []
            self.distance.append((event.xdata, event.ydata))
            print(self.keystrokes[event.key][1])
        if event.key == "D":
            if not self.distance:
                print(self.keystrokes[event.key][1])
                return
            else:
                self.distance.append((event.xdata, event.ydata))
                self.find_dist()
        if event.key == "G":
            #start blue polygon definition
            if not self.blue_polygon_active:
                self.polygon_vert_blue = []
                self.blue_polygon_active = True
            print(event.xdata, event.ydata)
            print(self.keystrokes[event.key][1])
            r = self._wcs.wcs_pix2world(event.xdata, event.ydata, 1)
            ra = float(r[0])
            dec = float(r[1])
            self.polygon_vert_blue.append([ra, dec])
            if len(self.polygon_vert_blue) > 2:
                self.show_lines([numpy.array(self.polygon_vert_blue[-2:]).T, ],
                                color='blue', linewidth=0.5)
            self.blue_polygon = True
            self.red_polygon = False
        if event.key == "g":
            #start red polygon definition
            if not self.red_polygon_active:
                self.polygon_vert_red = []
                self.red_polygon_active = True                
            print(event.xdata, event.ydata)
            print(self.keystrokes[event.key][1])
            r = self._wcs.wcs_pix2world(event.xdata, event.ydata, 1)
            ra = float(r[0])
            dec = float(r[1])
            self.polygon_vert_red.append([ra, dec])
            if len(self.polygon_vert_red) > 2:
                self.show_lines([numpy.array(self.polygon_vert_red[-2:]).T, ],
                                color='red', linewidth=0.5)
            self.blue_polygon = False
            self.red_polygon = True
        if event.key == "c":
            #close polygon definition
            if not self.polygon_vert_blue:
                if not self.polygon_vert_red:
                    print(self.keystrokes[event.key][1])
                    return
            if self.blue_polygon:
                self.polygon_vert_blue.append(self.polygon_vert_blue[0])
                self.pvert_blue = [numpy.array(self.polygon_vert_blue),]
                print(self.polygon_vert_blue)
                print(self.pvert_blue)
                #self.show_polygons(self.pvert_blue, edgecolor='b')
                self.blue_polygon_active = False
                self.polygon_average(self.pvert_blue, bluered='blue')
            else:
                self.polygon_vert_red.append(self.polygon_vert_red[0])
                self.pvert_red = [numpy.array(self.polygon_vert_red),]
                print(self.polygon_vert_red)
                print(self.pvert_red)
                #self.show_polygons(self.pvert_red, edgecolor='r')
                self.red_polygon_active = False
                self.polygon_average(self.pvert_red, bluered='red')
        if event.key == "h":
            #help text
            self.print_keystroke_help()

    def _get_hdu_title(self, num):
        if not self.extra_hdu_titles:
            title = 'HDU %d' % num
        else:
            title = self.extra_hdu_titles[num]
        return title

    def draw_spectrum(self, xdata, ydata):
        if not self.extra_hdus:
            print("Class was not instantiated with any extra HDU cubes")
            return
        self.specfigure.clf()
        xdata, ydata = map(int, (xdata, ydata))
        print(xdata, ydata)
        specax = self.specfigure.add_subplot(111)
        for i, hdu in enumerate(self.extra_hdus):
            spec = hdu.data[ydata, xdata, :]
            vel = getaxes(hdu.header, 1)
            title = self._get_hdu_title(i)
            specax.plot(vel, spec, drawstyle='steps-mid', label=title)
        specax.set_title('Spectrum at %.1f, %.1f' % (xdata, ydata))
        specax.legend(loc='best')
        self.specfigure.canvas.draw()    

    def draw_box_and_average(self):
        x = list(map(int, [f[0] for f in self.box]))
        y = list(map(int, [f[1] for f in self.box]))
        x.sort()
        y.sort()
        ax = self._figure.get_axes()[0]
        ax.hlines(y, x[0], x[1], colors='w')
        ax.vlines(x, y[0], y[1], colors='w')
        self.refresh()
        if self.avgfigure is None:
            self.avgfigure = mpl.figure()
        avgax = self.avgfigure.add_subplot(111)            
        for i, hdu in enumerate(self.extra_hdus):
            avgspec = numpy.nanmean(numpy.nanmean(hdu.data[y[0]:y[1], x[0]:x[1], :], axis=0), axis=0)
            vel = getaxes(hdu.header, 1)
            title = self._get_hdu_title(i)
            avgax.plot(vel, avgspec, drawstyle='steps-mid', label=title)
        avgax.set_title('Average of box (%.1f, %.1f) to (%.1f, %.1f)' % (x[0], y[0], x[1], y[1]))
        avgax.legend(loc='best')
        self.avgfigure.canvas.draw()

    def draw_box_and_sigma_average(self):
        x = list(map(int, [f[0] for f in self.box]))
        y = list(map(int, [f[1] for f in self.box]))
        x.sort()
        y.sort()
        ax = self._figure.get_axes()[0]
        ax.hlines(y, x[0], x[1], colors='w')
        ax.vlines(x, y[0], y[1], colors='w')
        self.refresh()
        avgfigure = mpl.figure()
        avgax = avgfigure.add_subplot(111)
        if self.blue_windows:
            b1, b2 = self.blue_windows
        else:
            raise SculptArgumentError('draw_box_and_sigma_average', 'You have to set blue windows to be able to do sigma averaged windows')
        if self.red_windows:
            r1, r2 = self.red_windows
        else:
            raise SculptArgumentError('draw_box_and_sigma_average', 'You have to set red windows to be able to do sigma averaged windows')
        
        for i, hdu in enumerate(self.extra_hdus):
            hdunew = cube_extract(hdu, [-1, -1, x[0], x[1], y[0], y[1]])
            vel = getaxes(hdu.header, 1)
            if self.blue:
                blue, brms = momentcube(hdunew, b1, b2, 
                                        returnrms=True, window=self.rms_window)
            else:
                red, rrms = momentcube(hdunew, r1, r2, 
                                       returnrms=True, window=self.rms_window)

            if self.blue:
                #blue-weighted
                yind, xind = numpy.where(blue.data > self.sigma_crit*brms.data/numpy.sqrt(2))
                avgspec = hdunew.data[yind, xind, :].mean(axis=0)
                ax.plot(xind+x[0], yind+y[0], 'bx')
                avgax.plot(vel, avgspec, drawstyle='steps-mid', label='blue wt %s' % self._get_hdu_title(i))
            else:
                #red-weighted
                yind, xind = numpy.where(red.data > self.sigma_crit*rrms.data/numpy.sqrt(2))
                avgspec = hdunew.data[yind, xind, :].mean(axis=0)
                ax.plot(xind+x[0], yind+y[0], 'rx')
                avgax.plot(vel, avgspec, drawstyle='steps-mid', label='red wt %s' % self._get_hdu_title(i))
        avgax.set_title("Average made from %s qualifying points" % len(xind))
        avgax.legend(loc='best')
        avgfigure.canvas.draw()
        self.refresh()

    def _get_points(self, hdu):
        self.ny, self.nx, nv = hdu.data.shape
        x = range(self.nx)
        y = range(self.ny)
        X, Y = numpy.meshgrid(x, y)
        x, y = X.flatten(), Y.flatten()
        xx, yy = xyad(hdu.header, x, y)
        return numpy.vstack((xx, yy)).T

    def polygon_average(self, pvert, bluered='blue'):
        """Averages within a polygon"""
        self.show_lines([numpy.array(pvert[0]).T,], edgecolor=bluered, linewidth=2)
        avgfigure = mpl.figure()
        avgax = avgfigure.add_subplot(111)
        for i, hdu in enumerate(self.extra_hdus):
            points = self._get_points(hdu)
            vel = getaxes(hdu.header, 1)
            #mask = points_inside_poly(points, pvert[0])
            ppath = path.Path(pvert[0])
            mask = ppath.contains_points(points) #inside_poly(points, pvert[0])
            navg = len(numpy.where(mask)[0])
            hdunew = hdu.data.copy()
            ny, nx, nv = hdunew.shape
            hdunew = hdunew.reshape((ny*nx), nv)
            poly_avg = hdunew[mask, :].mean(axis=0)
            avgax.plot(vel, poly_avg, drawstyle='steps-mid',
                       label='%s %s polygon' % (self._get_hdu_title(i), bluered))
        avgax.set_title("Average made from %s qualifying points" % navg)
        avgax.legend(loc='best')
        avgfigure.canvas.draw()
        self.refresh()

    def start_events(self):
        self.cid = self._figure.canvas.mpl_connect('key_press_event', self.on_keypress)
        

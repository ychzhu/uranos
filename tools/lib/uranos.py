# This is a copy of Corny/grains/uranos.py
# minor dependencies were commented out (3 comments) to run independently.

"""
CoRNy URANOS
    Functions specific for data processing of URANOS data
    based on: https://git.ufz.de/CRNS/cornish_pasdy/-/blob/master/corny/uranos.py
    version: 0.62
"""

import numpy as np
#from PIL import Image
from scipy.ndimage.measurements import label as scipy_img_label
import pandas
from glob import glob
import matplotlib.pyplot as plt

#from .corn import get_footprint, Wr, Wr_approx, sm2N_Koehli, sm2N, N2SM_Schmidt_single
#from .Schroen2017hess import get_footprint, Wr, Wr_approx
#from .Koehli2021fiw import sm2N_Koehli, sm2N, N2SM_Schmidt_single

class URANOS:

    variable_formats = dict({
        'Materials':'.0f',
        'Regions':'.0f',
        'region_id':'.0f',
        'center_mass':'',
        'center_geom':'',
        'area':'.1%',
        'SM':'.0%',
        'SM_diff':'+.0%',
        'Distance_min':'.1f',
        'Distance_com':'.1f',
        'Weights':'.1%',
        'Neutrons':'.0f',
        'Neutrons_diff':'+.0f',
        'Contributions':'.1%',
        'Contributions_diff':'+.1%',
        'Origins':'.1%',
        'Origins_err':'.1%',
        'Origins_below1':'.2%',
        'Density':'.2f',
        'freetext':''
    })
    
    variable_labels = dict({
        'Materials':      'Material code',
        'Regions':      'Regions (id)',
        'region_id':    'Regions (id)',
        'center_mass':  'Center of mass (grid)',
        'center_geom':  'Geometric center (grid)',
        'area':         'Area (%)',
        'SM':           'Soil Moisture (vol.%)',
        'SM_diff':           'Soil Moisture diff (vol.%)',
        'Distance_min': 'shortest Distance (m)',
        'Distance_com': 'Distance to center (m)',
        'Weights':      'Weights (%)',
        'Neutrons':     'estimated Neutrons (a.u.)',
        'Neutrons_diff':     'estimated Neutrons diff (a.u.)',
        'Contributions':'est. Contributions (%)',
        'Contributions_diff':'est. Contributions diff (%)',
        'Origins':      'sim. Contributions (%)',
        'Origins_err':      'sim. Contributions error (%)',
        'Origins_both':      'sim. Contributions (%)',
        'Density':      'sim. Neutron Density (a.u.)',
        'freetext':      'Materials'
    })
    
    def __init__(self, folder='', scaling=2, default_material=None, hum=5, press=1013, verbose=False):
        """
        Initialization
        """    
        self.verbose = verbose
        self.folder = folder
        self.scaling = scaling # one pixel in the data is x meters in reality
        self.center = (249, 249)
        self.default_material = default_material
        self.hum = hum
        self.press = press
        # Initial approximation with sm=20%, will be updated by materials2sm
        self.footprint_radius = 1#get_footprint(0.2, self.hum, self.press)
    
    #######
    # Input
    #######
    
    def read_materials(self, filename, scaling=None):
        """
        Read Material PNG image
        """
        I = Image.open(self.folder+filename)
        I = self.convert_rgb2grey(I)
        A = np.array(I)
        self.Materials = A
        self._idim = A.shape 
        self.center = ((self._idim[0]-1)/2, (self._idim[1]-1)/2)
        if not scaling is None:
            self.scaling = scaling
        print('Imported map `.Materials` (%d x %d), center at (%.1f, %.1f).' % (self._idim[0], self._idim[1], self.center[0], self.center[1]))
        print('  One pixel in the data equals %d meters in reality (%d x %d)' % (self.scaling, self._idim[0]*self.scaling, self._idim[1]*self.scaling))
        print('  Material codes: %s' % (', '.join(str(x) for x in np.unique(self.Materials))))
        if self.default_material is None:
            self.default_material = self.Materials[0,0]
            print('  Guessing default material: %d' % self.default_material)
        return(self)
        
    def read_origins(self, filepattern, pad=False):
        """
        Read URANOS origins matrix
        """
        U = None
        for filename in glob(self.folder + filepattern):
            print('  Reading %s' % filename)
            u = np.loadtxt(filename)
            if isinstance(U, np.ndarray):
                U += u
            else:
                U = u
        #U = np.loadtxt(self.folder+filename)
        if U is None:
            print('Error: no origin files found!')
            return(self)
        U = np.flipud(U)
        if pad:
            U = np.pad(U, ((0,1),(0,1)))
        self.Origins = U
        for i in self.region_data.index:
            self.region_data.loc[i, 'Origins'] = np.nansum(U[self.Regions==i]/np.nansum(U))
            self.region_data.loc[i, 'Origins_err'] = np.nansum(U[self.Regions==i]/np.nansum(U)) / np.sqrt(np.nansum(U[self.Regions==i]))

        print('Imported URANOS origins as `.Origins` (%d x %d).' % (U.shape[0], U.shape[1]))
        return(self)

    def read_density(self, filepattern, pad=False):
        """
        Read URANOS density matrix
        """
        U = None
        for filename in glob(self.folder + filepattern):
            print('  Reading %s' % filename)
            u = np.loadtxt(filename)
            if isinstance(U, np.ndarray):
                U += u
            else:
                U = u
            
        if U is None:
            print('Error: no files found!')
            return(self)
        
        U = np.flipud(U)
        if pad:
            U = np.pad(U, ((0,1),(0,1)))
        self.Density = U

        if hasattr(self, 'region_data'):
            for i in self.region_data.index:
                self.region_data.loc[i, 'Density'] = np.mean(U[self.Regions==i])
            self.region_data['Density'] /= self.region_data['Density'].max()

        print('Imported URANOS density map as `.Density` (%d x %d).' % (U.shape[0], U.shape[1]))
        return(self)

    def read_root(self, filepattern, show_vars=True):
        """
        Read URANOS root files
        """
        import uproot
        
        list_of_root_files = []
        for filename in glob(self.folder + filepattern):
            print('  Reading %s' % filename)
            rootfile = uproot.open(filename)
            if show_vars:
                print('Variables:', rootfile.keys())
            list_of_root_files.append(rootfile)

        if len(list_of_root_files)>1:
            self.Root = list_of_root_files
        elif len(list_of_root_files)==1:
            self.Root = list_of_root_files[0]
        else:
            print('Error: no root files found!')
            return(self)

        print('Imported URANOS root file as `.Root`. Consider reading vars with `read_root_var(var=...)`.')
        return(self)


    def read_root_var(self, var='detectorDistanceDepth2;1', alias=None, drop_zeros=True):
        """
        Convert a Root histogramm into a pandas DataFrame of coordinates
        Performance thanks to: Divakar <https://stackoverflow.com/a/41219731/2575273>
        """
        data = self.Root[var].to_numpy(flow=True)
        # convert histogramm to data frame
        m,n = data[0].shape
        R,C = np.mgrid[:m,:n]
        matrix = np.column_stack((C.ravel(),R.ravel(), data[0].ravel()))
        
        df = pandas.DataFrame(matrix, columns=['x','y','z'])
        df['x'] = df['x']*np.median(np.diff(data[2]))
        df['y'] = df['y']*np.median(np.diff(data[1]))
        
        if drop_zeros:
            df = df[df.z>0]

        if alias is None:
            alias = var

        if not hasattr(self, 'RootVar'):
            self.RootVar = dict()
        self.RootVar[alias] = df
        
        print('Imported URANOS root variable as `.RootVar[%s]` (%d points).' % (alias, len(df)))
        return(self)
        
    ############
    # Processing
    ############

    def D86_from_root(self, var='detectorDistanceDepth2;1', max_radius=25):
        """
        Calculate penetration depth D86 from root file
        """
        if not hasattr(self, 'RootVar') or not var in self.RootVar:
            self.read_root_var(var=var)

        df = self.RootVar[var]
        df = df[df.y < max_radius*1000]
        df = df.sort_values(['x'])
        df['z_cum'] = df.z.cumsum()
        
        D86 = df.loc[df.z_cum < 0.865*df.z_cum.max(),'x'].max()/1000
        return(D86)


    def material2sm(self, SM_gui=0.2):
        """
        Convert to soil moisture, assuming greyscale number between 2 and 139
        """
        
        self.SM = self.Materials/2/100
        self.SM[self.Materials == 100] = SM_gui  # water
        self.SM[self.Materials == 204] =  0.1  # concrete
        self.SM[self.Materials == 254] = 10.0  # water
        if self.Materials[(self.Materials<2)|(self.Materials>170)|(self.Materials!=204)|(self.Materials!=254)].any():
            print('Warning: some materials do not correspond to soil moisture, they are still naively converted using x/2/100.')
        print('Generated soil moisture map `.SM`, values: %s' % (', '.join(str(x) for x in np.unique(self.SM))))
        
        nearby_avg_sm = self.SM[self.m2grd(-25):self.m2grd(25), self.m2grd(-25):self.m2grd(25)]
        self.footprint_radius = get_footprint(nearby_avg_sm.mean(), self.hum, self.press)
        self.nearby_avg_sm = nearby_avg_sm
        print('Nearby avg. sm is %.2f +/- %.2f, updated footprint radius to %.0f m.' % (nearby_avg_sm.mean(), nearby_avg_sm.std(), self.footprint_radius))
        return(self)
    
    def generate_distance(self):
        """
        Distance matrix
        """
        D = np.zeros(shape=(self._idim[0],self._idim[1]))
        for i in range(self._idim[0]):
            for j in range(self._idim[1]):
                D[i,j] = np.sqrt(self.grd2m(i)**2 + self.grd2m(j)**2)
        self.Distance = D
        print('Generated distance map `.Distance`, reaching up to %.1f meters.' % (np.max(D)))
        return(self)
        
    
    def genereate_weights(self, approx=False, exclude_center=False):
        """
        Generate radial weights based on W_r
        """
        
        W = np.zeros(shape=(self._idim[0],self._idim[1]))
        for i in range(self._idim[0]):
            for j in range(self._idim[1]):
                r = 0.5 if self.Distance[i,j] == 0 else self.Distance[i,j] 
                if approx:
                    W[i,j] = Wr_approx(r)/r
                else:
                    W[i,j] = Wr(r, self.nearby_avg_sm.mean(), self.hum)/r
        if exclude_center:
            center_id = np.round(self.center).astype(int) 
            W[:,center_id[0]] = np.mean([W[:,center_id[0]-1], W[:,center_id[0]+1]])
            W[center_id[1],:] = np.mean([W[center_id[1]-1,:], W[center_id[0]+1,:]])
        W_sum = W.sum()
        self.Weights = W/W_sum
        print('Generated areal weights `.Weights`, ranging from %f to %f.' % (self.Weights.min(), self.Weights.max()))
        return(self)
        
    def find_regions(self, default_material=None):
        """
        Identidy connected regions based on the Materials map
        """
        if default_material is None:
            if self.default_material is None:
                default = np.unique(self.Materials)[0]
                print('Guessing default material code: %d' % default)
            else:
                default = self.default_material
        else:
            self.default_material = default_material
            default = default_material
        
        M = np.zeros(shape=(self._idim[0],self._idim[1]), dtype=np.uint8)
        for i in range(self._idim[0]):
            for j in range(self._idim[1]):
                M[i,j] = 0 if self.Materials[i,j] == default else 1
                
        #from scipy.ndimage.measurements import label
        L, ncomponents = scipy_img_label(M)
        self.Regions = L
        self.n_regions = ncomponents
        region_data = pandas.DataFrame(index=pandas.Series(np.arange(ncomponents+1), name='id'),
                                       columns=['Materials', 'center_mass', 'center_geom', 'area', 'SM',
                                                'Distance_min', 'Distance_com', 'Weights', 'Neutrons', 'Contributions', 'Origins', 'Density'])
        region_data['Regions'] = region_data.index
        for i in region_data.index:
            region_data.loc[i, 'Materials'] = np.median(self.Materials[self.Regions==i])
            region_data.at[i, 'center_mass'] = self._get_region_center(i, method='mass')
            region_data.at[i, 'center_geom'] = self._get_region_center(i, method='geom')
            region_data.loc[i, 'area'] = len(self.Regions[self.Regions==i]) / self._idim[0]/self._idim[1]
            if hasattr(self, 'SM'):
                region_data.loc[i, 'SM'] = np.median(self.SM[self.Regions==i])
            if hasattr(self, 'Weights'):
                region_data.loc[i, 'Weights'] = np.sum(self.Weights[self.Regions==i])
            if hasattr(self, 'Distance'):
                region_data.loc[i, 'Distance_min'] = np.min(self.Distance[self.Regions==i])
                region_data.loc[i, 'Distance_com'] = self.Distance[region_data.loc[i, 'center_mass'][0], region_data.loc[i, 'center_mass'][1]]
        
        self.region_data = region_data
        print('Found %d regions, mapped to `.Regions`, DataFrame generated as `.region_data`.' % self.n_regions)
        return(self)
    
    def estimate_neutrons(self, method='Koehli.2021', N0=1000, bd=1.43):
        """
        Estimate neutrons from the input soil moisture map
        """
        
        N = np.zeros(shape=(self._idim[0],self._idim[0]))
        for i in range(self._idim[0]):
            for j in range(self._idim[1]):
                if method=='Desilets.2010':
                    N[i,j] = sm2N(self.SM[i,j], N0, off=0.0, bd=bd)
                elif method=='Koehli.2021':
                    N[i,j] = sm2N_Koehli(self.SM[i,j], h=self.hum, off=0.0, bd=bd, func='vers2', method='Mar21_uranos_drf', bio=0)*N0/0.77
        
        self.Neutrons = N
        print('Estimated neutrons from soil moisture, `.Neutrons` (%.0f +/- %.0f)' % (N.mean(), N.std()))
        self.Contributions = self.Weights*N/np.sum(self.Weights*N)
        print('Estimated their signal contributions, `.Contributions`')
        
        for i in self.region_data.index:
            if not 'Neutrons_diff' in self.region_data.columns:
                self.region_data['Neutrons_diff'] = 0
            self.region_data.loc[i, 'Neutrons_diff'] = np.mean(N[self.Regions==i]) - self.region_data.loc[i, 'Neutrons']
            self.region_data.loc[i, 'Neutrons'] = np.mean(N[self.Regions==i])
            
            if not 'Contributions_diff' in self.region_data.columns:
                self.region_data['Contributions_diff'] = 0
            self.region_data.loc[i, 'Contributions_diff'] = np.sum(self.Contributions[self.Regions==i]) - self.region_data.loc[i, 'Contributions'] 
            self.region_data.loc[i, 'Contributions'] = np.sum(self.Contributions[self.Regions==i])
            
        return(self)

    def modify(self, Region=0, SM=None):
        if not SM is None:
            self.SM[self.Regions==Region] = float(SM)
            if not 'SM_diff' in self.region_data.columns:
                self.region_data['SM_diff'] = 0
            self.region_data.loc[self.region_data.Regions==Region, 'SM_diff'] = float(SM) - self.region_data.loc[self.region_data.Regions==Region, 'SM']
            self.region_data.loc[self.region_data.Regions==Region, 'SM'] = float(SM)
        return(self)
    
    ##################
    # Helper functions
    ##################
    
    def m2grd(self, m):
        g = m / self.scaling + self.center[0]
        if np.isscalar(g):
            return(int(np.round(g)))
        else:
            return(np.round(g).astype(int))
    def grd2m(self, g):
        m = g * self.scaling - self.center[0]*self.scaling
        if np.isscalar(m):
            return(int(np.round(m)))
        else:
            return(np.round(m).astype(int))
    
    def convert_rgb2grey(self, obj):
        """
        Convert an Array of Image from RGB to Greyscale
        """
        if isinstance(obj, np.ndarray):
            A = np.zeros(shape=(obj.shape[0],obj.shape[1]), dtype=np.uint8)
            for i in range(obj.shape[0]):
                for j in range(obj.shape[1]):
                    A[i,j] = np.mean(obj[i,j][0:3])
            return(A)

        elif isinstance(obj, Image.Image):
            img = obj.convert('L')
            return(img)
        
        else:
            print('Error: provide either a numpy array of RGB arrays or an RGB image.')
            
    def _get_region_center(self, region_id, method='mass'):
        """
        Get the center of a region
        """
        indices = np.where(self.Regions==region_id)
        if method=='geom':
            yrange = (indices[0].min(), indices[0].max())
            xrange = (indices[1].min(), indices[1].max())
            center = (0.5*xrange[0]+0.5*xrange[1], 0.5*yrange[0]+0.5*yrange[1])
        elif method=='mass':
            center = (indices[1].mean(), indices[0].mean())
        return([int(np.round(x)) for x in center])
    
    ########
    # Output
    ########
    
    def plot(self, ax=None, image='SM', annotate=None, overlay=None, fontsize=10, title=None, contour=False,
             regions=None, extent=500, cmap='Greys', vmax_factor=2, x_marker=None, cross_alpha=0.5,
             label_offset=(0,0)):
        """
        Plot map, annotate, and overlay
        Avoid annotation with annotate='none'
        """
        
        # If no regions are provided, show all regions
        if regions is None and hasattr(self, 'region_data'):
            regions = np.arange(len(self.region_data))
        
        if annotate is None:
            annotate = image
            
        try:
            Var = getattr(self, image)
        except:
            print('Error: %s is no valid attribute.' % image)
        
        if ax is None:
            fig, ax = plt.subplots(figsize=(5,5))
        
        if title is None:
            title = 'Map of %s' % self.variable_labels[image]
            if annotate != image:    title += '\nAnnotation: %s' % self.variable_labels[annotate]
            if overlay=='Origins':   title += '\nOverlay: sim. Neutron Origins (x)'
            if not x_marker is None: title += '\n\n'
        else:
            title = title
        ax.set_title(title, fontsize=fontsize)
        ax.imshow(Var, interpolation='none', cmap=cmap, vmax=np.max(Var)*vmax_factor)
        if contour:
            ax.contour(self.Regions, levels=0, colors='k', linewidths=1, antialiased=True)

        lox, loy = label_offset
        if regions:
            for i in regions:
                mask = (self.Regions == i)
                dataset = self.region_data.loc[i]
                
                cmx, cmy = dataset['center_mass']
                coords = (cmx+lox*i, cmy+loy*i)
                
                # all this clutter only for having 2 digits for origins < 1...
                fmt = self.variable_formats[annotate] if annotate in self.variable_formats else '%s'
                
                if annotate in ['Origins','Origins_err','Origins_both']:
                    if dataset['Origins'] < 0.01:
                        fmt     = self.variable_formats['Origins_below1']
                        fmt_err = self.variable_formats['Origins_below1']
                    else:
                        fmt     = self.variable_formats['Origins']
                        fmt_err = self.variable_formats['Origins_err']
                elif annotate == 'Contributions':
                    if dataset['Contributions'] < 0.01:
                        fmt = self.variable_formats['Origins_below1']
                    
                if annotate=='Origins_both':
                    ax.annotate('{0:{1}}\n$\pm${2:{3}}'.format(dataset['Origins'], fmt, dataset['Origins_err'], fmt_err),
                        coords, fontsize=fontsize, ha='center', va='center')
                else:
                    ax.annotate('{0:{1}}'.format(dataset[annotate], fmt),
                        coords, fontsize=fontsize, ha='center', va='center')

        #plt.title(r'$\theta_1=5\,\%$, $\theta_2=10\,\%$, $R=200\,$m')
        # Tick format in meters
        grid100 = np.arange(-extent*0.8, extent*0.8+1, extent*0.2)
        ax.set_xticks(self.m2grd(grid100), ['%.0f' % x for x in grid100])
        ax.set_yticks(self.m2grd(grid100), ['%.0f' % x for x in grid100[::-1]])
        # Zoom to extend
        ax.set_xlim(self.m2grd(grid100-extent*0.2).min(), self.m2grd(grid100+extent*0.2).max())
        ax.set_ylim(self.m2grd(grid100+extent*0.2).max(), self.m2grd(grid100-extent*0.2).min())
        # Tick markers inside
        ax.tick_params(axis="y", direction='in', length=4)
        ax.tick_params(axis="x", direction='in', length=4)
        # Ticks an all sides
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        # Central lines
        #ax.axhline(m2grd(0), color='black', lw=1, ls='--', alpha=0.3)
        #ax.axvline(m2grd(0), color='black', lw=1, ls='--', alpha=0.3)
        ax.plot(self.m2grd(0),self.m2grd(0), color='black', marker='+', ms=10, mew=1, alpha=cross_alpha)

        # Distance marker
        if not x_marker is None:
            ax.annotate('$R=%.0f\,$m' % x_marker, xy=(self.m2grd(x_marker), 0), xytext=(self.m2grd(x_marker), -35),
                        ha='center', va='center', arrowprops={'arrowstyle': 'wedge'}) 
            # for arrows, see https://matplotlib.org/stable/gallery/text_labels_and_annotations/fancyarrow_demo.html

        if not overlay is None:
            try:
                Var = getattr(self, overlay)
            except:
                print('Error: %s is no valid attribute.' % overlay)
            x,y = np.meshgrid(np.arange(Var.shape[1]), np.arange(Var.shape[0]))
            Var[np.where(Var < 1)] = np.nan
            #ax.imshow(X.Origins, interpolation='none')
            ax.scatter(x,y, c=Var, s=24, cmap='autumn', alpha=0.3, marker='x')
        
        # Footprint circle
        fpc = plt.Circle((self.m2grd(0), self.m2grd(0)), self.footprint_radius/self.scaling, fc='none', ec='black', ls='--', alpha=0.3)
        ax.add_patch(fpc)
        
        return(ax)


    def histogram(self, ax=None, var='SM'):
        """
        Plot a histogram of a certain variable in .region_data 
        """
        if not var in self.region_data:
            print('Error: Variable is not in .region_data')
            return()
        
        if ax is None:
            fig, ax = plt.subplots(figsize=(5,3))
            
        self.region_data[var].plot.barh(color='black', ax=ax)
        ax.grid(axis='x', color='black', alpha=0.3)
        ax.invert_yaxis()
        ax.set_yticks(range(len(self.region_data)+1))
        ax.set_ylabel('Region number')
        ax.set_xlabel(self.variable_labels[var])
        return(ax)
    
    def average_sm(self, N0=1000):
        """
        Calculate different approaches to get average CRNS soil moisture
        """
        print('Average soil moisture seen by the CRNS detector:')
        print('{0:.1%}             field mean (naive approach)'.format(self.SM.mean()))
        print('{0:.1%} SM-weighted field mean (lazy approach)'.format((self.Weights*self.SM).sum()))
        print('{0:.1%}  N-weighted field mean (correct approach)'.format(N2SM_Schmidt_single((self.Weights*self.Neutrons).sum()/N0*0.77, bd=1.43, hum=self.hum)))

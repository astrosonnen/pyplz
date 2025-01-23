import numpy as np
from imageSim import SBModels


def read_config(filename):

    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    config = {'data_dir':'./', 'mask_dir': None, 'output_dir':'./', 'sps_model_dir': None, 'filters': None, \
              'zeropoints': None, \
              'filename': None, 'filter_prefix': '', 'filter_suffix': '', 'science_tag':'_sci.fits', 'err_tag':'_var.fits', 'err_type': 'VAR', 'psf_tag':'_psf.fits', \
              'rmax': None, 'Nsteps':300, 'Nwalkers':30, 'burnin':None, 'maskname':None, \
              'rgbscheme': 'LINEAR', 'rgbcuts': None, 'rgbscales': None, 'outname': None, 'modeltype': 'standard', 'reference_band': None}

    preamble = True

    i = 0
    while preamble and i < len(lines):
        if '#' in lines[i] and 'MODELS' in lines[i]:
            preamble = False
        else:
            line = lines[i].split('#')[0].split()
            if len(line) > 0:
                parname = line[0].split(':')[0]
                if parname in config:
                    parval = lines[i].split('#')[0].split(':')[1].split('\n')[0].strip()
                    config[parname] = parval

        i += 1

    filtlist = []
    filternames = config['filters'].split(',')
    for name in filternames:
        filtlist.append(name.strip())
    config['filters'] = filtlist

    config['rgbscheme'] = config['rgbscheme'].upper()

    if config['rgbscheme'] == 'LINEAR':
        config['rgbscales'] = None
        if config['rgbcuts'] is not None:
            cutlist = []
            cuts = config['rgbcuts'].split('')
            for cut in cuts:
                cutlist.append(float(cut))

            config['rgbcuts'] = cutlist
        else:
            config['rgbcuts'] = (99., 99., 99.)

    elif config['rgbscheme'] == 'M16':
        config['rgbcuts'] = None
        scalelist = []
        if config['rgbscales'] is not None:
            scales = config['rgbscales'].split(',')
            for scale in scales:
                scalelist.append(float(scale))
        else:
            for name in filternames:
                scalelist.append(1.)
        config['rgbscales'] = scalelist
    else:
        raise NameError("rgbscheme can only be 'LINEAR' or 'M16'")

    config['zeropoints'] = np.array(config['zeropoints'].split(','), dtype='float')
    config['Nsteps'] = int(config['Nsteps'])
    config['Nwalkers'] = int(config['Nwalkers'])
    if config['burnin'] is not None:
        config['burnin'] = int(config['burnin'])
    config['rmax'] = float(config['rmax'])

    light_components = []
    lens_components = []
    source_components = []

    while i < len(lines):

        line = lines[i].split()
        if len(line) > 0:
            if line[0] == 'light_model':
                model_class = line[1]

                if model_class in SBModels.parlists:
                    parnames = []
                    for parname in SBModels.parlists[model_class]:
                        parnames.append(parname)
                else:
                    sys.exit("model class '%s' is not supported"%model_class)

                comp = {'class':model_class, 'pars':{}}

                if config['modeltype'] == 'SED':
                    sed_class = line[2]
                    comp['sed'] = sed_class
                    sps_model = None
                    if sed_class == 'colorpars':
                        config['colors'] = []
                        for band in config['filters']:
                            if band != config['reference_band']:
                                config['colors'].append('%s-%s'%(band, config['reference_band']))
                        parnames += config['colors']
                    elif sed_class == 'template':
                        parnames += ['zs', 'tn']
                    elif sed_class == 'sps':
                        sps_model = line[3]
                        parnames += ['redshift', 'age', 'tau', 'logZ', 'logtau_V']
                    else:
                        df

                npars = len(parnames)

                foundpars = 0
                j = 1

                while foundpars < npars and j+i < len(lines):
                    line = lines[j+i].split()
                    if lines[j+i][0] != '#' and len(line) > 0:
                        if line[0] in parnames:
                            foundpars += 1
                            par = line[0]
                            link = None
                            if len(line) > 6:
                                link = line[6]
                            tmp_par = {'value': float(line[1]), 'low': float(line[2]), 'up': float(line[3]), \
                           'step': float(line[4]), 'var': int(line[5]), 'link':link}
                            comp['pars'][par] = tmp_par
                    j += 1

                i += j

                if foundpars < npars:
                    print('not all parameters found!')
                else:
                    light_components.append(comp)

            elif line[0] == 'source_model':
                model_class = line[1].strip()

                comp = {'class':model_class, 'pars':{}}

                if model_class in ['Sersic']:
                    parnames = []
                    for parname in SBModels.parlists[model_class]:
                        parnames.append(parname)
                else:
                    df

                if config['modeltype'] == 'SED':
                    sed_class = line[2]
                    comp['sed'] = sed_class
                    sps_model = None
                    if sed_class == 'colorpars':
                        config['colors'] = []
                        for band in config['filters']:
                            if band != config['reference_band']:
                                config['colors'].append('%s-%s'%(band, config['reference_band']))
                        parnames += config['colors']
                    elif sed_class == 'template':
                        parnames += ['zs', 'tn']
                    elif sed_class == 'sps':
                        sps_model = line[3]
                        parnames += ['redshift', 'age', 'tau', 'logZ', 'logtau_V']
                    else:
                        df

                npars = len(parnames)

                foundpars = 0
                j = 1

                while foundpars < npars and j+i < len(lines):
                    line = lines[j+i].split()
                    if lines[j+i][0] != '#' and len(line) > 0:
                        if line[0] in parnames:
                            foundpars += 1
                            par = line[0]
                            link = None
                            if len(line) > 6:
                                link = line[6]
                            tmp_par = {'value': float(line[1]), 'low': float(line[2]), 'up': float(line[3]), \
                           'step': float(line[4]), 'var': int(line[5]), 'link':link}
                            comp['pars'][par] = tmp_par
                    j += 1

                i += j

                if foundpars < npars:
                    print('not all parameters found!')
                else:
                    source_components.append(comp)

            elif 'lens_model' in line[0]:
                model_class = line[1].strip()

                if model_class == 'Powerlaw':
                    npars = 6
                    parnames = ['x', 'y', 'q', 'pa', 'b', 'eta']
                else:
                    df

                comp = {'class': model_class, 'pars':{}}

                foundpars = 0
                j = 1

                while foundpars < npars and j+i < len(lines):
                    line = lines[j+i].split()
                    if lines[j+i][0] != '#' and len(line) > 0:
                        if line[0] in parnames:
                            foundpars += 1
                            par = line[0]
                            link = None
                            if len(line) > 6:
                                link = line[6]
                            tmp_par = {'value': float(line[1]), 'low': float(line[2]), 'up': float(line[3]), \
                           'step': float(line[4]), 'var': int(line[5]), 'link':link}
                            comp['pars'][par] = tmp_par
                    j += 1

                i += j

                if foundpars < npars:
                    print('not all parameters found!')
                else:
                    lens_components.append(comp)

            else:
                i += 1
        else:
            i += 1

    config['light_components'] = light_components
    config['source_components'] = source_components
    config['lens_components'] = lens_components

    return config


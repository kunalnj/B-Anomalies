import uproot
import numpy as np
import pandas as pd
from numpy.linalg import solve

def DOCA_point(df,mu_plus_point,mu_minus_point,bs_point,bs_dir,mu_plus_dir,mu_minus_dir,n_plus,n_minus):
    '''
    Get point on Bs LOF where DOCA occurs
    '''
    coords_plus=np.zeros([len(mu_plus_point[:,0]),3])

    for i in range(len(mu_plus_point[:,0])):

        rhs=mu_plus_point[i,:]-bs_point[i,:]
        lhs=np.array([bs_dir[i,:], -mu_plus_dir[i,:], n_plus[i,:]]).T
        X=solve(lhs,rhs)

        coords_plus[i,:]=bs_point[i,:]+X[0]*bs_dir[i,:] # point on Bs LoF where we have the CDA. Basically the reconstructed decay vertex of the tau

    df['tauplus_closest_x'] = coords_plus[:,0]
    df['tauplus_closest_y'] = coords_plus[:,1]
    df['tauplus_closest_z'] = coords_plus[:,2]

    coords_minus=np.zeros([len(mu_minus_point[:,0]),3])

    for i in range(len(mu_minus_point[:,0])):

        rhs=mu_minus_point[i,:]-bs_point[i,:]
        lhs=np.array([bs_dir[i,:],-mu_minus_dir[i,:],n_minus[i,:]]).T
        X=solve(lhs,rhs)

        coords_minus[i,:]=bs_point[i,:]+X[0]*bs_dir[i,:] # point on Bs LoF where we have the CDA. Basically the reconstructed decay vertex of the tau

    df['tauminus_closest_x'] = coords_minus[:,0]
    df['tauminus_closest_y'] = coords_minus[:,1]
    df['tauminus_closest_z'] = coords_minus[:,2]

def FD(df):
    '''
    Analytical calculation of FD
    '''
    tauplus_closest_vars = ['tauminus_closest_x', 'tauminus_closest_y', 'tauminus_closest_z']
    tauminus_closest_vars = ['tauminus_closest_x', 'tauminus_closest_y', 'tauminus_closest_z']

    tau_start = df[['phi3_VX', 'phi3_VY', 'phi3_VZ']].to_numpy()
    tauplus_end = df[tauplus_closest_vars].to_numpy() # THE CLOSEST DISTANCE POINT
    tauminus_end = df[tauminus_closest_vars].to_numpy() # THE CLOSEST DISTANCE POINT

    # tauplus_end = df[['tauplus_VX', 'tauplus_VY', 'tauplus_VZ']].to_numpy() # THE CLOSEST DISTANCE POINT
    # tauminus_end = df[['tauminus_VX', 'tauminus_VY', 'tauminus_VZ']].to_numpy() # THE CLOSEST DISTANCE POINT

    diff_plus = tauplus_end - tau_start
    diff_minus = tauminus_end - tau_start

    sign_plus = +1 # inner1d(bs_dir,diff_plus)/abs(inner1d(bs_dir,diff_plus))
    sign_minus = +1 # inner1d(bs_dir,diff_minus)/abs(inner1d(bs_dir,diff_minus))

    mag_plus = np.sum(np.abs(diff_plus)**2,axis=-1)**(1./2)
    mag_minus = np.sum(np.abs(diff_minus)**2,axis=-1)**(1./2)

    df['tauplus_analytic_FD'] = mag_plus * sign_plus # added the sign!!! Important if want to knwo which way the tau went basically
    df['tauminus_analytic_FD'] = mag_minus * sign_minus # added the sign!!! Important if want to knwo which way the tau went basically

def IP(df,bs_point,mu_plus_point,mu_minus_point,mu_plus_dir,mu_minus_dir):
    '''
    Analytic calc of IP
    '''
    n_mu_plus = mu_plus_dir/np.linalg.norm(mu_plus_dir, axis = 1, keepdims = True)
    n_mu_minus = mu_minus_dir/np.linalg.norm(mu_minus_dir, axis = 1, keepdims = True)

    d_along_mu_plus = np.abs(((mu_plus_point - bs_point)*n_mu_plus).sum(1))
    d_along_mu_minus = np.abs(((mu_minus_point - bs_point)*n_mu_minus).sum(1))

    diff_plus = np.linalg.norm( (mu_plus_point - bs_point) , axis = 1)
    diff_minus = np.linalg.norm( (mu_minus_point - bs_point) , axis = 1)

    IP_plus = np.sqrt(np.abs(diff_plus**2 - d_along_mu_plus**2))
    IP_minus = np.sqrt(np.abs(diff_minus**2 - d_along_mu_minus**2))

    df['IP_analytic_plus'] = IP_plus
    df['IP_analytic_minus'] = IP_minus

def root_to_df(fname,drop_na = True, analytic_doca = True):
    '''
    Create DataFrame from root file, append IP, DOCA, FD
    '''

    if '/' and 'Bs2KKmumu' in fname:
        file = uproot.open('{}'.format(fname))
        decay_tuple = file['DecayTree']
        col_names = decay_tuple.keys()
        columns = [col.decode("utf-8") for col in col_names]
        # df = pd.DataFrame(decay_tuple.arrays(columns, flatten=None))
        # df.columns = columns
        #df = decay_tuple.pandas.df(columns, flatten=True)
        Bs = decay_tuple.pandas.df("Bs*", flatten=True)
        phi = decay_tuple.pandas.df("phi*", flatten=True)
        muminus = decay_tuple.pandas.df("mu_minus*", flatten=True)
        muplus = decay_tuple.pandas.df("mu_plus*", flatten=True)
        tauminus = decay_tuple.pandas.df("tauminus*", flatten=True)
        tauplus = decay_tuple.pandas.df("tauplus*", flatten=True)
        Kminus = decay_tuple.pandas.df("K_minus*", flatten=True)
        Kplus = decay_tuple.pandas.df("K_plus*", flatten=True)
        # dfs = [Bs, phi, dist, mu, Kplus, Kminus]
        df = Bs.join(phi, how="outer").join(muplus, how="outer").join(muminus, how="outer").join(tauplus, how="outer").join(tauminus, how="outer").join(Kminus, how="outer").join(Kplus, how="outer")
        if drop_na:
            df.dropna(inplace=True)
        '''
        Analytic calculation of Distance of Closest Approach (DOCA)
        '''
        bs_point_1 = df[['Bs_OWNPV_X', 'Bs_OWNPV_Y', 'Bs_OWNPV_Z']].to_numpy() #on Bs LOF
        bs_point_2 = df[['Bs_ENDVERTEX_X', 'Bs_ENDVERTEX_Y', 'Bs_ENDVERTEX_Z']].to_numpy() #on Bs LOF
        mu_plus_point = df[['mu_plus_REFPX', 'mu_plus_REFPY', 'mu_plus_REFPZ']].to_numpy() #on mu + LOF
        #mu_plus_point = df[['mu_plus_OWNPV_X','mu_plus_OWNPV_Y','mu_plus_OWNPV_Z']].applymap(lambda x: x[0]).to_numpy() #on mu + LOF
        mu_minus_point = df[['mu_minus_REFPX', 'mu_minus_REFPY', 'mu_minus_REFPZ']].to_numpy() #on mu - LOF
        #mu_minus_point = df[['mu_minus_OWNPV_X','mu_minus_OWNPV_Y','mu_minus_OWNPV_Z']].applymap(lambda x: x[0]).to_numpy() #on mu - LOF

        bs_dir = bs_point_1 - bs_point_2
        mu_plus_dir = df[['mu_plus_AtVtx_PX', 'mu_plus_AtVtx_PY', 'mu_plus_AtVtx_PZ']].to_numpy()
        mu_minus_dir = df[['mu_minus_AtVtx_PX', 'mu_minus_AtVtx_PY', 'mu_minus_AtVtx_PZ']].to_numpy()

        n_plus = np.cross(bs_dir,mu_plus_dir)
        n_plus = n_plus/np.linalg.norm(n_plus, axis = 1, keepdims = True)
        n_minus = np.cross(bs_dir,mu_minus_dir)
        n_minus = n_minus/np.linalg.norm(n_minus, axis = 1, keepdims = True)

        DOCA_plus = np.abs(((mu_plus_point - bs_point_1)*n_plus).sum(1))
        DOCA_minus = np.abs(((mu_minus_point - bs_point_1)*n_minus).sum(1))

        df['DOCA_analytic_plus'] = DOCA_plus
        df['DOCA_analytic_minus'] = DOCA_minus

        np.savetxt('columns/{}_columns.txt'.format(fname.split('.')[0].split('/')[-1]),df.columns.to_numpy(),fmt='%s')

    else:
        file = uproot.open('~/code/root_files/{}'.format(fname))
    
        # decay_tuple = file['DecayTuple/DecayTuple']
        # col_names = decay_tuple.keys()
        # columns = [col.decode("utf-8") for col in col_names]
        # df = pd.DataFrame(decay_tuple.arrays(columns))
        # df.columns = columns

        decay_tuple = file['DecayTuple/DecayTuple']
        col_names = decay_tuple.keys()
        columns = [col.decode("utf-8") for col in col_names]
        # df = pd.DataFrame(decay_tuple.arrays(columns, flatten=None))
        # df.columns = columns
        #df = decay_tuple.pandas.df(columns, flatten=True)
        Bs = decay_tuple.pandas.df("Bs*", flatten=True)
        phi = decay_tuple.pandas.df("phi3*", flatten=True)
        dist = decay_tuple.pandas.df("DOCA*", flatten=True)
        muminus = decay_tuple.pandas.df("muminus*", flatten=True)
        muplus = decay_tuple.pandas.df("muplus*", flatten=True)
        tauminus = decay_tuple.pandas.df("tauminus*", flatten=True)
        tauplus = decay_tuple.pandas.df("tauplus*", flatten=True)
        Kminus = decay_tuple.pandas.df("Kminus*", flatten=True)
        Kplus = decay_tuple.pandas.df("Kplus*", flatten=True)
        # dfs = [Bs, phi, dist, mu, Kplus, Kminus]
        df = Bs.join(phi, how="outer").join(dist, how="outer").join(muplus, how="outer").join(muminus, how="outer").join(tauplus, how="outer").join(tauminus, how="outer").join(Kminus, how="outer").join(Kplus, how="outer")

        if drop_na: df.dropna(inplace=True)

        df['DOCA_mu_plus_err'] = df['DOCA_mu_plus']/np.sqrt(df['DOCA_chi2_mu_plus'])
        df['DOCA_mu_minus_err'] = df['DOCA_mu_minus']/np.sqrt(df['DOCA_chi2_mu_minus'])

        '''
        Analytic calculation of Distance of Closest Approach (DOCA)
        '''
        bs_point = df[['phi3_VX','phi3_VY','phi3_VZ']].to_numpy() #on Bs LOF
        mu_plus_point = df[['muplus_RefPoint_X','muplus_RefPoint_Y','muplus_RefPoint_Z']].to_numpy() #on mu + LOF
        mu_minus_point = df[['muminus_RefPoint_X','muminus_RefPoint_Y','muminus_RefPoint_Z']].to_numpy() #on mu - LOF

        bs_dir = df[['phi3_PX', 'phi3_PY','phi3_PZ']].to_numpy()
        mu_plus_dir = df[['muplus_PX','muplus_PY','muplus_PZ']].to_numpy()
        mu_minus_dir = df[['muminus_PX','muminus_PY','muminus_PZ']].to_numpy()

        n_plus = np.cross(bs_dir,mu_plus_dir)
        n_plus = n_plus/np.linalg.norm(n_plus, axis = 1, keepdims = True)
        n_minus = np.cross(bs_dir,mu_minus_dir)
        n_minus = n_minus/np.linalg.norm(n_minus, axis = 1, keepdims = True)

        DOCA_plus = np.abs(((mu_plus_point - bs_point)*n_plus).sum(1))
        DOCA_minus = np.abs(((mu_minus_point - bs_point)*n_minus).sum(1))

        df['DOCA_analytic_plus'] = DOCA_plus
        df['DOCA_analytic_minus'] = DOCA_minus

        DOCA_point(df,mu_plus_point,mu_minus_point,bs_point,bs_dir,mu_plus_dir,mu_minus_dir,n_plus,n_minus)
        FD(df)
        IP(df,bs_point,mu_plus_point,mu_minus_point,mu_plus_dir,mu_minus_dir)
        
        np.savetxt('columns/{}_columns.txt'.format(fname.split('.')[0]),df.columns.to_numpy(),fmt='%s')

    return df



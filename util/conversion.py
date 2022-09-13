import sys, glob, uuid
import numpy as np
import h5py as h5
import ROOT as rt
import subprocess as sub
from util.config import GetNPars, GetJetConfig, GetInvisiblesFlag, GetSignalFlag
from util.calcs import DeltaR2, EPxPyPzToPtEtaPhiM, EPzToRap, PtEtaPhiMToPxPyPzE, PtEtaPhiMToEPxPyPz, AdjustPhi, EPxPyPzToM
from util.fastjet import BuildFastjet
from util.qol_utils.qol_util import printProgressBarColor, RN
from util.qol_utils.pdg import pdg_plotcodes, pdg_names, FillPdgHist
from util.conv_utils.utils import ExtractHepMCParticles, InitFastJet, PrepDataBuffer, ExtractHepMCEvents,PrepDelphesArrays, PrepH5File, PrepIndexRanges
from util.particle_selection import IsNeutrino
# --- FASTJET IMPORT ---
# TODO: Can this be done more nicely?
fastjet_dir = BuildFastjet(j=8)
fastjet_dir = glob.glob('{}/**/site-packages'.format(fastjet_dir),recursive=True)[0]
if(fastjet_dir not in sys.path): sys.path.append(fastjet_dir)
import fastjet as fj
# ----------------------

# Convert HepMC / Delphes-ROOT files into our final output -- possibly performing jet clustering along the way.
class Processor:
    def __init__(self, use_delphes=False):
        self.delphes = use_delphes

        self.prefix_level1 = '\tClustering jets & preparing data:'
        self.suffix = 'Complete'
        self.bl = 50

        self.jet_config, self.jetdef = InitFastJet()

        self.outdir = ''

        self.hist_filename = 'hists.root'
        self.diagnostic_plots = True
        self.InitializeHists()

    def SetOutputDirectory(self,outdir):
        self.outdir = outdir

    def SetHistFilename(self,name):
        self.hist_filename = name

    def SetDiagnosticPlots(self,flag):
        self.diagnostic_plots = flag

    def InitializeHists(self):
        self.hists = []
        self.pdg_hist = rt.TH1D(RN(),'Selected final-state particles;Particle;Count',47,0.,47.)
        self.hists.append(self.pdg_hist)
        self.jet_truth_dr_hists = {} # Entries get initialized and filled in FillDataBuffer

        self.jet_hists = { # Jet kinematic histograms. Initialized/filled in FillDataBuffer.
            'pt'  : rt.TH1D(RN(), 'Jet p_{T};p_{T} [GeV];Count', 500,0.,1000.),
            'eta' : rt.TH1D(RN(), 'Jet #eta;#eta;Count',         200,-4.,4.),
            'phi' : rt.TH1D(RN(), 'Jet #phi;#phi;Count',         200,-np.pi,np.pi),
            'm'   : rt.TH1D(RN(), 'Jet m;m [GeV];Count',         250,0.,500.),
            'e'   : rt.TH1D(RN(), 'Jet E;E [GeV];Count',         500,0.,1000.),
            'et'  : rt.TH1D(RN(), 'Jet E_{T};E [GeV];Count',     500,0.,1000.)
        }

        self.pmu_sum_hists = {} # Ditto, for the sum of jet constituents.
        for key,hist in self.jet_hists.items():
            h = rt.TH1D(hist)
            h.SetName(RN())
            h.SetTitle(h.GetTitle().replace('Jet','#sum constituents'))
            self.pmu_sum_hists[key] = h

        self.jet_diff_hists = {
            'dm': rt.TH1D(RN(),  ';(m_{consts} - m_{j}) / m_{j};Count'      , 400,-1.,1.),
            'de': rt.TH1D(RN(),  ';(E_{consts} - E_{j}) / E_{j};Count'      , 400,-1.,1.),
            'dpt': rt.TH1D(RN(), ';(p_{consts} - p_{T,j}) / p_{T,j};Count', 400,-1.,1.),
            'det': rt.TH1D(RN(), ';(E_{consts} - E_{T,j}) / E_{T,j};Count', 400,-1.,1.),
        } # Kinematic differences between the jets, and the sum of their constituents. (Sanity check!)

        self.jet_diff_hists_2d = {
            'dm': rt.TH2D(RN(),  ';m_{j};(m_{consts} - m_{j}) / m_{j};Count',         250,0.,500. ,400,-1.,1.),
            'de': rt.TH2D(RN(),  ';E_{j};(E_{consts} - E_{j}) / E_{j};Count',         500,0.,1000.,400,-1.,1.),
            'dpt': rt.TH2D(RN(), ';p_{T,j};(p_{consts} - p_{T,j}) / p_{T,j};Count', 500,0.,1000.,400,-1.,1.),
            'det': rt.TH2D(RN(), ';E_{T,j};(E_{T,consts} - E_{T,j}) / E_{T,j};Count', 500,0.,1000.,400,-1.,1.),
        } # Kinematic differences between the jets, and the sum of their constituents. (Sanity check!)

    def OutputHistograms(self):
        for key,h in self.jet_truth_dr_hists.items():
            self.hists.append(h)

        for key,h in self.pmu_sum_hists.items():
            self.hists.append(h)

        for key,h in self.jet_hists.items():
            self.hists.append(h)

        for key,h in self.jet_diff_hists.items():
            self.hists.append(h)

        for key,h in self.jet_diff_hists_2d.items():
            self.hists.append(h)

        rt.gStyle.SetOptStat(0)
        offset = 10 # TODO: This is a bit hacky
        hist_filename = '{}/{}'.format(self.outdir,self.hist_filename)
        f = rt.TFile(hist_filename,'UPDATE')

        # We separately handle the pdg_hist (first in list), since it requires
        # modifying the axis labels after drawing.
        for i,h in enumerate(self.hists):
            canvas_name = 'c_{}'.format(offset + i)
            hist_name = 'h_{}'.format(offset + i)
            c = rt.TCanvas(canvas_name,'',800,600)

            if(i == 0): # pdg_hist

                pdg_names_inv = {}
                for key,val in pdg_plotcodes.items():
                    code = val
                    name = pdg_names[key]
                    pdg_names_inv[code] = name

                n = h.GetNbinsX()
                xaxis = h.GetXaxis()
                for j in range(n):
                    if(j == 0): name = ''
                    else:
                        try: name = pdg_names_inv[j]
                        except: name = ''
                    xaxis.SetBinLabel(j+1,name)
                xaxis.SetLabelSize(2.0e-2)
                rt.gPad.SetGrid(1,0)
            else: rt.gPad.SetGrid()

            if(type(h) in [rt.TH1F, rt.TH1D, rt.TH1I, rt.TH1S]): h.Draw('HIST')
            else:
                h.Draw('COLZ')
                rt.gPad.SetRightMargin(0.2)
            rt.gPad.SetLogy()
            c.Draw()
            f.cd() # unnecessary
            # c.Write(canvas_name)
            h.Write(hist_name)
        f.Close()
        return

    def Process(self, final_state_files, truth_files=None, h5_file=None, nentries_per_chunk=1e4, verbosity=0, separate_truth_particles=True):
        if(type(final_state_files) == list): final_state_files = ['{}/{}'.format(self.outdir,x) for x in final_state_files]
        else: final_state_files = '{}/{}'.format(self.outdir,final_state_files)

        if(h5_file is None):
            if(type(final_state_files) == list): h5_file = final_state_files[0]
            else: h5_file = final_state_files
            if(self.delphes): h5_file =  h5_file.replace('*','').replace('.root','.h5')
            else: h5_file =  h5_file.replace('*','').replace('.root','.h5')
        else:
            h5_file = '{}/{}'.format(self.outdir,h5_file)
        npars = GetNPars()
        n_truth = npars['n_truth']

            # Note: It's important that the Delphes files and truth HepMC files line up!
        #       The way we usually do this, it will be guaranteed.
        if(type(final_state_files) == str): final_state_files = glob.glob(final_state_files,recursive=True)
        if(truth_files is None):
            if(self.delphes): truth_files = [x.replace('.root','_truth.hepmc') for x in final_state_files]
            else: truth_files = [x.replace('.hepmc','_truth.hepmc') for x in final_state_files]
        else:
            truth_files = ['{}/{}'.format(self.outdir,x) for x in truth_files]

        if(self.delphes):
            delphes_arr,var_map = PrepDelphesArrays(final_state_files)
            nentries = len(delphes_arr)
        else:
            ## Extract final-state truth particle info from the HepMC files.
            final_state_events, nentries = ExtractHepMCEvents(final_state_files,get_nevents=True)

        # Also extract truth particle info from the HepMC files.
        truth_events = ExtractHepMCEvents(truth_files)

        nentries_per_chunk = int(nentries_per_chunk)
        data = PrepDataBuffer(nentries_per_chunk,separate_truth_particles)

        # Prepare the HDF5 file.
        dsets = PrepH5File(h5_file,nentries,data)

        # Some indexing preparation for writing in chunks.
        start_idxs,stop_idxs,ranges = PrepIndexRanges(nentries,nentries_per_chunk)
        nchunks = len(ranges)

        if(verbosity == 1): printProgressBarColor(0,nchunks, prefix=self.prefix_level1, suffix=self.suffix, length=self.bl)

        for i in range(nchunks):
            # Clear the buffer (for safety).
            for key in data.keys(): data[key][:] = 0

            if(self.delphes):
                # Get the momenta of the final-state particles, first in cylindrical coordinates.
                # We will be setting m=0 for all these particles.
                momenta = {}
                types = var_map.keys() # TODO: Can this be placed elsewhere? Renamed?
                for delphes_type in types:
                    momenta[delphes_type] = {
                        x:delphes_arr[var_map[delphes_type][x]][start_idxs[i]:stop_idxs[i]]
                        for x in ['pt','eta','phi']
                    }
            else:
                ## Get final-state particles.
                npar_fs = 1e4 # TODO: This is some hardcoded max number of final-state particles. Should be plenty.
                final_state_particles = ExtractHepMCParticles(final_state_events[start_idxs[i]:stop_idxs[i]],npar_fs)

            # Get the truth particles.
            truth_particles = ExtractHepMCParticles(truth_events[start_idxs[i]:stop_idxs[i]],n_truth)

            prefix_nzero = int(np.ceil(np.log10(nchunks))) + 1
            prefix_level2 = '\tClustering jets & preparing data for chunk {}/{}:'.format(str(i+1).zfill(prefix_nzero),nchunks)
            if(verbosity==2): printProgressBarColor(0,ranges[i], prefix=prefix_level2, suffix=self.suffix, length=self.bl)

            # For each event we must combine tracks and neutral hadrons, perform jet clustering on them,
            # select a single jet (based on user criteria), and select that jet's leading constituents.
            for j in range(ranges[i]):
                data['is_signal'][j] = GetSignalFlag()

                if(self.delphes):
                    pt  = np.concatenate([momenta[x]['pt'][j].to_numpy() for x in types])
                    eta = np.concatenate([momenta[x]['eta'][j].to_numpy() for x in types])
                    phi = np.concatenate([momenta[x]['phi'][j].to_numpy() for x in types])
                    m   = np.zeros(pt.shape)

                    vecs = PtEtaPhiMToPxPyPzE(pt,eta,phi,m) # convert 4-vectors to (px, py, pz, e) for jet clustering.
                    px = vecs[:,0]
                    py = vecs[:,1]
                    pz = vecs[:,2]
                    e  = vecs[:,3]
                    pdg = None # No meaning to PDG codes if we're looking at detector-level reconstruction

                else:
                    l = len(final_state_particles[j])

                    # Optionally exclude invisibles (neutrinos). If using Delphes, this should be handled by Delphes internally.
                    if(GetInvisiblesFlag()):
                        visibles = np.arange(l)
                    else:
                        invisibles = np.array([IsNeutrino(final_state_particles[j][k]) for k in range(l)])
                        visibles = np.where(invisibles == 0)[0]
                        del invisibles

                    px = np.array([final_state_particles[j][k].momentum.px for k in visibles])
                    py = np.array([final_state_particles[j][k].momentum.py for k in visibles])
                    pz = np.array([final_state_particles[j][k].momentum.pz for k in visibles])
                    e  = np.array([final_state_particles[j][k].momentum.e  for k in visibles])
                    pdg  = np.array([final_state_particles[j][k].pid  for k in visibles])

                self.SelectFinalStateParticles(px,py,pz,e, pdg, self.jetdef, truth_particles, data, j, separate_truth_particles)

                if(verbosity==2): printProgressBarColor(j+1,ranges[i], prefix=prefix_level2, suffix=self.suffix, length=self.bl)

            # We have now filled a chunk, time to write it.
            with h5.File(h5_file, 'a') as f:
                for key in dsets.keys():
                    dset = f[key]
                    dset[start_idxs[i]:stop_idxs[i]] = data[key][:ranges[i]]
            if(verbosity == 1): printProgressBarColor(i+1,nchunks, prefix=self.prefix_level1, suffix=self.suffix, length=self.bl)

        if(self.diagnostic_plots): self.OutputHistograms()
        return h5_file

    # --- Utility functions for the class ---
    def InitFastJet(self):
        fj.ClusterSequence.print_banner() # Get the Fastjet banner out of the way.
        jet_config = GetJetConfig()
        jetdef = fj.JetDefinition(fj.antikt_algorithm, jet_config['jet_radius'])
        return jet_config,jetdef

    def ClusterJets(self,vecs, jetdef, jet_config):
        pj = [fj.PseudoJet(*x) for x in vecs]
        jets = jetdef(pj)

        # Apply optional minimum jet pT cut.
        jet_pt = np.array([jet.pt() for jet in jets])
        jet_indices = np.linspace(0,len(jets)-1,len(jets),dtype=np.dtype('i8'))[jet_pt >= jet_config['jet_min_pt']]
        jets = [jets[i] for i in jet_indices]

        # Apply optional maximum |eta| cut.
        jet_eta = np.array([jet.eta() for jet in jets])
        jet_indices = np.linspace(0,len(jets)-1,len(jets),dtype=np.dtype('i8'))[np.abs(jet_eta) <= jet_config['jet_max_eta']]
        jets = [jets[i] for i in jet_indices]
        return jets

    def FetchJetConstituents(self,jet,n_constituents):
        # pt,eta,phi,m = np.hsplit(np.array([[x.pt(), x.eta(), x.phi(), x.m()] for x in jet.constituents()]),4)
        pt,eta,phi,m,px,py,pz,e = np.hsplit(np.array([[x.pt(), x.eta(), x.phi(), x.m(), x.px(),x.py(),x.pz(),x.e()] for x in jet.constituents()]),8)
        pt = pt.flatten() # use this for sorting
        eta = eta.flatten()
        phi = phi.flatten()
        m = m.flatten()
        px = px.flatten()
        py = py.flatten()
        pz = pz.flatten()
        e  =  e.flatten()

        # Sort by decreasing pt, and only keep leading constituents.
        sorting = np.argsort(-pt)
        l = int(np.minimum(n_constituents,len(pt)))

        pt = pt[sorting][:l]
        # eta = eta[sorting][:l]
        # phi = phi[sorting][:l]
        # m  =  m[sorting][:l]
        px = px[sorting][:l]
        py = py[sorting][:l]
        pz = pz[sorting][:l]
        e  =  e[sorting][:l]

        vecs = np.array([e,px,py,pz],dtype=np.dtype('f8')).T

        # vecs = PtEtaPhiMToEPxPyPz(pt=pt, # Vector will be of format (E, px, py, pz)
        #                             eta=eta,
        #                             phi=phi,
        #                             m=m
        #                         )
        return vecs

    def SelectFinalStateParticles(self,px,py,pz,e, pdg, jetdef, truth_particles, data, j, separate_truth_particles):
        jet_config = GetJetConfig()
        n_constituents = GetNPars()['jet_n_par']
        jet_sel = jet_config['jet_selection']

        if(jet_sel is None):
            vecs = np.vstack((e,px,py,pz)).T # note energy is given as 0th component
            self.FillDataBuffer(data,j,vecs,None,truth_particles,separate_truth_particles)
        else:
            vecs = np.vstack((px,py,pz,e)).T # note energy is given as last component, needed for fastjet
            jets = self.ClusterJets(vecs,jetdef,jet_config)

            njets = len(jets)
            if(njets == 0):
                data['is_signal'][j] = -1 # Using -1 to mark as "no event". (To be discarded.)
                return

            # selected_jet_idx = jet_config['jet_selection'](truth=truth_particles[j], jets=jets, use_hepmc=True)
            selected_jet_idx = jet_sel(truth=truth_particles[j], jets=jets)
            if(selected_jet_idx < 0):
                data['is_signal'][j] = -1 # Using -1 to mark as "no event". (To be discarded.)
                return
            jet = jets[selected_jet_idx]
            # Get the constituents of our selected jet.
            vecs = self.FetchJetConstituents(jet,n_constituents)
            self.FillDataBuffer(data,j,vecs,jet,truth_particles,separate_truth_particles)

        # If PDG code information was provided, let's determine the codes of the selected particles.
        # TODO: This will involve looping -> will be intensive. Is there a better way?
        if(pdg is not None and self.diagnostic_plots):
            pdg_matching_tolerance = 1.0e-4
            selected_pdgs = np.zeros(pdg.shape,dtype=bool)
            candidate_particles = np.vstack((e,px,py,pz)).T
            for i,particle in enumerate(vecs):
                for j,candidate_particle in enumerate(candidate_particles):
                    if(selected_pdgs[j] == 1): continue
                    match_value = np.dot(candidate_particle-particle,candidate_particle-particle)
                    if match_value < pdg_matching_tolerance:
                        selected_pdgs[j] = 1
                        break
            selected_pdgs = pdg[selected_pdgs]
            FillPdgHist(self.pdg_hist,selected_pdgs)
        return

    def FillDataBuffer(self,data_buffer,j,vecs,jet,truth_particles,separate_truth_particles=False):
        npars = GetNPars()
        n_truth = npars['n_truth']
        l = vecs.shape[0]
        data_buffer['Nobj'][j] = l
        data_buffer['Pmu'][j,:l,:] = vecs
        data_buffer['truth_Nobj'][j] = len(truth_particles[j])
        data_buffer['truth_Pdg'][j,:] = [par.pid for par in truth_particles[j]]
        data_buffer['truth_Pmu'][j,:,:] = [
            [par.momentum.e, par.momentum.px, par.momentum.py, par.momentum.pz]
            for par in truth_particles[j]
        ]

        if(jet is None):
            data_buffer['jet_Pmu'][j,:] = 0.
            data_buffer['jet_Pmu_cyl'][j,:] = 0.
        else:
            data_buffer['jet_Pmu'][j,:]     = [jet.e(), jet.px(), jet.py(), jet.pz()]
            data_buffer['jet_Pmu_cyl'][j,:] = [jet.pt(), jet.eta(), AdjustPhi(jet.phi()), jet.m()]

            if(self.diagnostic_plots):
                # Histogram the dR between the jet and each of the truth-level particles.
                for i, pid in enumerate(data_buffer['truth_Pdg'][j,:]):
                    if(pid not in self.jet_truth_dr_hists.keys()):
                        self.jet_truth_dr_hists[pid] = rt.TH1D('jet_truth_dr_{}'.format(pid),'',200,0.,2.)
                        title = '#Delta R #left( jet, {}_{} #right)'.format(pdg_names[pid],'{truth}')
                        self.jet_truth_dr_hists[pid].SetTitle(title)
                        self.jet_truth_dr_hists[pid].GetXaxis().SetTitle('#Delta R')
                        self.jet_truth_dr_hists[pid].GetYaxis().SetTitle('Count')
                    h = self.jet_truth_dr_hists[pid]
                    dr = np.sqrt(DeltaR2(truth_particles[j][i].momentum.eta(),truth_particles[j][i].momentum.phi(), data_buffer['jet_Pmu_cyl'][j,1],data_buffer['jet_Pmu_cyl'][j,2]))
                    h.Fill(dr)

        if(separate_truth_particles):
            for k in range(n_truth):
                data_buffer['truth_Pmu_{}'.format(k)][j,:] = [truth_particles[j][k].momentum.e,truth_particles[j][k].momentum.px,truth_particles[j][k].momentum.py,truth_particles[j][k].momentum.pz]

        # Let's make some more diagnostic plots, to look at the jet kinematics, and the kinematics of the sum of jet constituents.
        if(self.diagnostic_plots):

            # Kinematic histograms for sum of Pmu (jet constituents).
            s = np.sum(vecs,axis=0) # E, px, py, pz
            s_cyl = EPxPyPzToPtEtaPhiM(*s)
            s_pt,s_eta,s_phi,s_m = s_cyl
            s_e = s[0]
            y = EPzToRap(s[0],s[3])
            s_et = s[0]/np.cosh(y)
            self.pmu_sum_hists['pt'].Fill(s_pt)
            self.pmu_sum_hists['eta'].Fill(s_eta)
            self.pmu_sum_hists['phi'].Fill(s_phi)
            self.pmu_sum_hists['m'].Fill(s_m)
            self.pmu_sum_hists['e'].Fill(s_e)
            self.pmu_sum_hists['et'].Fill(s_et)

            if(jet is not None):

                # Kinematic histograms for jets.
                jet_pt = data_buffer['jet_Pmu_cyl'][j,0]
                jet_eta = data_buffer['jet_Pmu_cyl'][j,1]
                jet_phi = data_buffer['jet_Pmu_cyl'][j,2]
                jet_m = data_buffer['jet_Pmu_cyl'][j,3]
                jet_e = jet.e()
                jet_et = jet.Et()

                self.jet_hists['pt'].Fill(jet_pt)
                self.jet_hists['eta'].Fill(jet_eta)
                self.jet_hists['phi'].Fill(jet_phi)
                self.jet_hists['m'].Fill(jet_m)
                self.jet_hists['e'].Fill(jet_e)
                self.jet_hists['et'].Fill(jet_et)

                # Kinematic difference histograms.
                delta_e = (s_e - jet_e) / jet_e
                delta_m = (s_m - jet_m) / jet_m
                delta_pt = (s_pt - jet_pt) / jet_pt
                delta_et = (s_et - jet_et) / jet_et

                self.jet_diff_hists['de'].Fill(delta_e)
                self.jet_diff_hists['dm'].Fill(delta_m)
                self.jet_diff_hists['dpt'].Fill(delta_pt)
                self.jet_diff_hists['det'].Fill(delta_et)

                self.jet_diff_hists_2d['de'].Fill(jet_e,delta_e)
                self.jet_diff_hists_2d['dm'].Fill(jet_m,delta_m)
                self.jet_diff_hists_2d['dpt'].Fill(jet_pt,delta_pt)
                self.jet_diff_hists_2d['det'].Fill(jet_et,delta_et)

        return

# Remove any "failed events". They will be marked with a negative signal flag.
def RemoveFailedFromHDF5(h5_file, cwd=None):
    if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)

    fname_tmp = h5_file.replace('.h5','_{}.h5'.format(str(uuid.uuid4())))

    f = h5.File(h5_file,'r')
    keys = list(f.keys())
    shapes = [f[key].shape for key in keys]
    N = shapes[0][0] # number of events

    keep_indices = f['is_signal'][:] >= 0
    if(np.sum(keep_indices) == N): return

    g = h5.File(fname_tmp,'w')

    for i,key in enumerate(keys):
        g.create_dataset(key, data = f[key][:][keep_indices],compression='gzip')

    f.close()
    g.close()

    sub.check_call(['rm',h5_file])
    sub.check_call(['mv',fname_tmp, h5_file])
    return

# Split an HDF5 file into training, testing and validation samples.
def SplitHDF5(h5_file, split_ratio = (7,2,1), train_name=None, val_name=None, test_name=None, cwd=None):
    if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)
    file_dir = '/'.join(h5_file.split('/')[:-1])

    if(train_name is None):
        train_name = 'train.h5'
        if(file_dir != ''): train_name = '{}/{}'.format(file_dir,train_name)
    if(val_name is None):
        val_name   =   'valid.h5'
        if(file_dir != ''): val_name = '{}/{}'.format(file_dir,val_name)
    if(test_name is None):
        test_name  =  'test.h5'
        if(file_dir != ''): test_name = '{}/{}'.format(file_dir,test_name)

    f = h5.File(h5_file,'r')
    keys = list(f.keys())
    shapes = [f[key].shape for key in keys]
    N = shapes[0][0] # number of events

    split_ratio = np.array(split_ratio)
    tot = np.sum(np.array(split_ratio))
    split_ratio = split_ratio / tot
    split_ratio[-1] = 1. - np.sum(split_ratio[:-1])

    names = [train_name, val_name, test_name]
    n = np.array(N * split_ratio,dtype=int)
    diff = np.sum(n) - N
    n[-1] -= diff

    rng = np.random.default_rng()
    index_list = np.array(range(N),dtype=int)
    rng.shuffle(index_list,axis=0)
    indices = []
    start = 0
    stop = 0
    for i in range(len(names)):
        stop += n[i]
        idxs = index_list[start:stop]
        start = stop
        indices.append(idxs)

    for i in range(len(names)):
        g = h5.File(names[i],'w')
        for j,key in enumerate(keys):
                # Putting the data from f into memory (numpy array)
                # since h5py doesn't like non-ordered indices.
                g.create_dataset(key, data=f[key][:][indices[i]],compression='gzip')
        g.close()
    f.close()
    return
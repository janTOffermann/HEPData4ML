import sys, os, glob
import  ROOT as rt
import numpy as np
import h5py as h5
import argparse as ap

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-i','--input',type=str,required=True,help='Input file.')
    parser.add_argument('-r','--reference',type=str,required=True,help='Reference HDF5 file, input will be rescaled to match jet pT distribution of these.')
    parser.add_argument('-o','--output',type=str,required=True,help='Output file.')
    parser.add_argument('-pt1','--pt_min',type=float,help='Minimum pt for rescaling.',default=200.)
    parser.add_argument('-pt2','--pt_max',type=float,help='Maximum pt for rescaling.',default=1000.)
    parser.add_argument('-nbins','--nbins',type=int,help='Number of pt bins.',default=100)
    args = vars(parser.parse_args())

    infile = args['input']
    reference = args['reference']
    outfile = args['output']
    pt_min = args['pt_min']
    pt_max = args['pt_max']
    nbins = args['nbins']

    outdir = '/'.join(outfile.split('/')[:-1])
    os.makedirs(outdir,exist_ok=True)

    binning = np.linspace(pt_min,pt_max,nbins+1)

    f = h5.File(infile,'r')
    input_pt = f['jet_Pmu_cyl'][:,0]

    g = h5.File(reference,'r')
    ref_pt = g['jet_Pmu_cyl'][:,0]

    g.close()

    h_ref_pt = rt.TH1D('ref_pt','',len(binning)-1,binning)
    h_inp_pt = rt.TH1D('inp_pt','',len(binning)-1,binning)
    for entry in ref_pt:
        h_ref_pt.Fill(entry)
    for entry in input_pt:
        h_inp_pt.Fill(entry)

    h_scale = h_ref_pt / h_inp_pt
    nbins = h_scale.GetNbinsX()

    for i in range(nbins):
        if(h_scale.GetBinContent(i+1) > 1):
            h_scale.SetBinContent(i+1,1.)

    h_inp_pt_scaled = h_scale * h_inp_pt

    print("Original  # of events: ",h_inp_pt.Integral())
    print("Rescaled  # of events: ",h_inp_pt_scaled.Integral())
    print("Reference # of events: ",h_ref_pt.Integral())
    print("\tratio = {:.2e}".format(h_inp_pt_scaled.Integral() / h_ref_pt.Integral()))

    c = rt.TCanvas('c','c',800,600)
    h_ref_pt.Draw('HIST')
    h_ref_pt.SetLineColor(rt.kRed)
    h_ref_pt.SetLineWidth(2)
    h_inp_pt_scaled.Draw('HIST SAME')

    pt_bin_number = np.array([h_inp_pt.FindBin(x) for x in input_pt],dtype=int)
    # print(pt_bin_number)
    rng = np.random.default_rng()
    randoms = np.array([rng.uniform() for x in pt_bin_number],dtype=float)
    # scale_factors = np.array([h_scale.GetBinContent(int(bn)) for bn in pt_bin_number])
    flag = np.where(np.array([randoms[i] <= h_scale.GetBinContent(int(bn)) for i,bn in enumerate(pt_bin_number)],dtype=bool))[0]

    keys = list(f.keys())
    meta_keys = list(f.attrs.keys())
    g = h5.File(outfile,'w')

    for key in keys:
        g.create_dataset(key,data=f[key][:][flag],compression='gzip',compression_opts=9)

    for key in meta_keys:
        g.attrs[key] = f.attrs[key]

    f.close()
    g.close()

    # h_scale.Draw('HIST')
    # h_ref_pt.Draw('HIST')
    # h_inp_pt.Draw('HIST SAME')
    c.Draw()
    c.SaveAs("{}/pt.pdf".format(outdir))

if(__name__ == '__main__'):
    main(sys.argv)
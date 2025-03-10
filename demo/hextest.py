
import os
import platform
import pathlib

def check_and_make_dir(din):
    d = din.resolve()
    dnpy = (din / "*.npy").resolve() 
    if not os.path.exists(d):
        os.makedirs(d)
    else:
        if platform.system()=='Windows':
            print(f"del {dnpy}")
            os.system(f"del {dnpy}")
        else:
            print(f"rm {dnpy}")
            os.system(f"rm {dnpy}")

def check_for_output( d, f):
    check = os.path.exists( (d / f).resolve() )
    if check:
        r = "Passed"
    else:
        r = "Failed"
    return r

def print_result( sname, result):
    print("\n\n*******************************")
    print(f"****** {sname} test : {result}" )
    print("*******************************\n")

test_diff = True
test_corr = True
test_mask = True
test_padf = True
test_plot = True
test_dc   = False

npatterns = 8   #change this to 1000 to reproduce the paper results; set it to a low number like 6 to quickly check the code is working 

outdir = pathlib.Path( "./output")

# diffraction calculation
if test_diff:
    check_and_make_dir( outdir / "diff")
    print("\n Performing diffract.py test")
    os.system(f"python ../diffract.py --config ./configs/config_hex_diff.txt --npatterns {npatterns}")
    test = check_for_output( outdir / "diff", "hex_2D_0.npy")
    #print("\n diffract.py test :", test)
    print_result( "diffract.py", test)

# correlation calculation
if test_corr:
    check_and_make_dir( outdir / "corr")
    print("\n Performing difftocorr.py test")
    os.system(f"python ../difftocorr.py --config ./configs/config_hex_corr.txt --npatterns {npatterns} --outputsum True")
    test = check_for_output( outdir / "corr", "hex_a_correlation_sum.npy")
    print("\n difftocorr.py test :", test)
    print_result( "difftocorr.py", test)

# correlation calculation
if test_dc:
    check_and_make_dir( outdir / "dc")
    print("\n Performing diffract_and_correlate.py test")
    os.system("python ../diffract_and_correlate.py --config ./configs/config_hex_dc.txt")
    test = check_for_output( outdir / "dc", "hexdc_a_correlation_sum.npy")
    print("\n diffract_and_correlate.py test :", test)
    print_result( "diffract_and_correlate.py", test)

# maskcorr calculation
if test_mask:
    print("\n Performing maskcorr.py test")
    os.system("python ../maskcorr.py --config ./configs/config_hex_mask.txt")
    test = check_for_output( outdir / "corr", "hex_ab_correlation_sum_sintheta.npy")
    print("\n maskcorr.py test passed")
    print_result( "maskcorr.py", test)

# padf calculation
if test_padf:
    check_and_make_dir( outdir / "padf")
    print("\n Performing corrtopadf.py test")
    os.system("python ../corrtopadf.py --config ./configs/config_hex_padf.txt") 
    test = check_for_output( outdir / "padf", "hex_ab_padf.npy")
    print("\n corrtopadf.py test :", test)
    print_result( "corrtopadf.py", test)

# plot calculation
if test_plot:
    check_and_make_dir( outdir / "figs")
    print("\n Performing plotfxs3d.py test")
    os.system("python ../plotfxs3d.py --config ./configs/config_hex_plot.txt")
    os.system("python ../plotfxs3d.py --config ./configs/config_hex_plot.txt --convolve True --suffix convolved")
    os.system("python ../plotfxs3d.py --config ./configs/config_hex_corrplot.txt")
    os.system("python ../plotfxs3d.py --config ./configs/config_hex_corrplot.txt --fname ./output/corr/hex_ab_correlation_sum_sintheta.npy")
    test = check_for_output( outdir / "figs", "hex_ab_padf_test_reqr.npy")
    print("\n plotfxs3d.py test :", test)
    print_result( "plotfxs3d.py", test)

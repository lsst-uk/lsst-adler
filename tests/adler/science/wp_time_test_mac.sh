time astscript-radial-profile -q -h1 /Users/jrobinson/lsst-adler/tests/data/ztf_Didymos-system-barycenter20065803_20221201402616_000567_zr_c10_o_q1_scimrefdiffimg.fits.fz -a 0.0,36.0 --measure=sum,mean,median,sigclip-mean,sigclip-std --center=32.0,32.5 -o ./wedge_out_0.txt -R 32.0; 
time cat ./wedge_out_0.txt; 
time rm ./wedge_out_0.txt;

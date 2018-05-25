
addprocs(Sys.CPU_CORES-nprocs());
nprocs()

push!(LOAD_PATH,"library");
using myFun
using NK_tools
# using PyPlot

# Change the following three lines:
xdim,ydim,zdim,tdim = 128,128,20,28;
filenameReal = string("data/real_x128_y128_z20_d28.float32")
filenameImag = string("data/imag_x128_y128_z20_d28.float32")

# Loading data
fid = open(filenameReal,"r");
dataReal = read(fid, Float32,xdim,ydim,zdim,tdim);
close(fid);
fid = open(filenameImag,"r");
dataImag = read(fid, Float32,xdim,ydim,zdim,tdim);
close(fid);
data = dataReal + complex(0,1)*dataImag;

# Phase correction
hann(x) =  flipdim(0.5 * (1 - cos.(2π*(1:x)/x)),1); 
hann2d = zeros(xdim,ydim);
hann2d[div(xdim,2)-div(xdim,4)+1:div(xdim,2)+div(xdim,4),div(ydim,2)-div(ydim,4)+1:div(ydim,2)+div(ydim,4)] = hann(div(xdim,2))*hann(div(ydim,2))';
data_pc = zeros(Complex64,size(data));
for cnt = 1:tdim
    tmp1 = data[:,:,:,cnt];
    tmp10 = copy(tmp1);
    k1 = qift(tmp1);
    k1s = k1.*hann2d;
    tmp1s = qft(k1s);
    tmp2 = tmp10.*exp.(-1im*angle.(tmp1s));
    data_pc[:,:,:,cnt] = tmp2;
end

# Denoising
@time data_sm = denoisingDWI_complex(data_pc);
# figure(1,figsize=(6,2));imal(flipdim(abs.(data[:,:,10,12:14]),1),3,1);
# figure(2,figsize=(6,2));imal(flipdim(abs.(data_sm[:,:,10,12:14]),1),3,1);

# Saving the results
# change the filename below:
filename = string("data/","magnitude_x128_y128_z20_d28_denoised.float32");
fid = open(filename,"w");
write(fid, data_sm);
close(fid);

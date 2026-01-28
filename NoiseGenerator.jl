using FFTW

function NoiseGenerator(N, alpha)
    # length of time series, parameter of power law (pink noise: alpha=1)
    # Fourier frequencies
    f = range(0, stop=pi, length=N÷2+1)[2:end-1]
    # power law
    ff = 1.0 ./ f .^ alpha
    # real of Fourier transform
    RW = sqrt.(0.5 .* ff) .* randn(N÷2-1)
    # imaginary of Fourier transform
    IW = sqrt.(0.5 .* ff) .* randn(N÷2-1)
    fR = complex.([randn(); RW; randn(); reverse(RW)], [0; IW; 0; -1 .* reverse(IW)])

    # complex numbers to be retransformed
    # sequence of frequencies: 0,2pi/N, 2*2pi/N,...,pi,...,2pi-1/N
    # two frequencies ff1 = pi + c , ff2 = pi - c shall be complex conjugates
    # frequencies at 0 or pi shall have imaginary=0

    # transform to times
    reihe = fft(fR, 1) / N
    # series as reals
    return(real(reihe))
end

# Example usage:
 result = NoiseGenerator(1000,-2.0)
 plot(result,label="α=3.0")

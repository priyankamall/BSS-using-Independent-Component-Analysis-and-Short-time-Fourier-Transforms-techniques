# BSS-using-Independent-Component-Analysis-and-Short-time-Fourier-Transforms-techniques
This experiment focuses on the task of separating a mixture of independent (non-gaussian)
sound signals using the technique of independent component analysis (ICA). The idea of 
blind source separation (BSS) is to determine the original signals given the linear mixtures.
In ICA, one assumes that the mixtures and the sources are random variables (with zero mean) 
instead of proper time signals. This allows the final mixed signals to be considered as 
samples of these random variables. To denote these mathematically, consider s to be the 
matrix of source signals (where each row consists of the samples of the source); let the 
matrix A be the mixing matrix, which produces the mixed signals x by applying a linear 
transformation on the sources. We have that: 
x = As (1) 
Now, the objective of BSS is to recover the source signals s. given only the matrix x.
More explicitly, using ICA, the objective is to find a matrix W such that:           
s = Wx (2) 
Note that although the goal is to recover the source signals as accurately as possible, 
it is not necessary for the sources (rows) to be ordered in the same manner as in the
original matrix s in Equation 2. In other words, the first row of the recovered matrix s′ 
may not be the same as that of the signal in the first row of s, but it has
to be very close to one of the rows of s. 


When separating mixtures of instantaneously mixed signals, independent component analysis
works very well, but this ideal situation is rarely found in the real-world. In practical
applications, the environment distorts audio signals by adding echoes, reflections, and 
ambient noise. Additionally independent component analysis, in its purest form, assumes 
that source signals do not have any propagation delay, which is an assumption that cannot 
be applied in this case. Recording sources from two microphones placed in different locations 
will inevitably introduce propagation delays, so the blind source separation method used 
must also consider this issue. To solve problems detailed above, the blind source separation 
problem will be redefined in the time-frequency domain. By taking the short-time Fourier
transform (STFT) of the audio inputs, we can represent the inputs as the following.



![Your Repository's Stats](https://github-readme-stats.vercel.app/api/top-langs/?username=Your_GitHub_Username&theme=blue-green)

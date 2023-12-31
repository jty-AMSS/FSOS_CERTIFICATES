# FSOS_CERTIFICATES

This is the a tool to compute the FSOS cetrificates of function on finite abelian groups, which is based on SOS relaxation.  

## Dependencies
- Matlab 
- [SDPNAL+](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/)

We have test it in Matlab 2016b and Matlab 2020b
## Usage
The Numerical experiments of this paper can be performed by "Numerical_experiments.m"

### construct a function on fintie abelian group  
$$f:G \mapsto \mathbb{Z},~G=\prod_{i=1}^{k}  C_{n_i},$$
with 

$$f=\sum_{j=1}^{s}a_jx^{\alpha_j}, ~x_i \in  C_{n_i}$$
Let n be a 1-dim array such that n(i)=n_i, then set
```
f=CZ(n);
for j=1:s
f(alpha_j)=a_j;
end
```
### compute polynomial/rational FSOS
Let $f$ be a function on fintie abelian group with range in set S, one can use 
```
 [Q,lb,Index,A,B]=short_proof_poly_fsos(f,d,S,ADMmaxiterCoeff);
```
and
```
 [U,V,S,T]=short_proof_rational_fsos(f,d,S,ADMmaxiterCoeff)
```
to compute the polynomial/rational FSOS, and use
```
 [err,H]=CheckFSOS(f,Index,Q)
```
to check the polynomial FSOS cetrificate ,where $H$ is a cell of functions on group $G$,  err equals to $\sum_{i} |H_i|^2$.
```
[err,p,q,H_1,H_2]=CheckFSOS_rational(f,S,T,U,V) 
```
to check the rational FSOS, where err equals to $p-q \cdot f$, and  $H_1$, $H_2$  are cells of functions on group $G$, with  $p=\sum_{i}|H_{1,i}|^2$, $q=\sum_{i}|H_{2,i}|^2$. 

# Demos
DEMO.m are the computation of examples shown in Section 3-Section 5 

BGP.m are the examples in Section 6.1

Demo_3sat.m are the examples in Section 6.2.1

GraphExample.m and GraphExample_SAT.m are the examples in Section 6.3.1

PHP_Example.m are the examples in Section 6.3.2

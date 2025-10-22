from estimator import *

print("For one-more-MISIS scheme:\n")

power = 25		# bit size of the modulus
q = 2**power
deg = 768		# degree of the ring R_q
k = deg*2        	# degree * dim
bkz_d = 2*k     	# BKZ dim is degree * dim * 2
w = 16			# weight of sparse LWE error
gamma = 2		# rejection parameter

n_sigma_0 = 1.5*sqrt(q) *2*sqrt(deg)	# norm of sigma_0
s = gamma * n_sigma_0 * (1+sqrt(3*w))	# standard deviation of signature
bound = 2*s*sqrt(deg)			# norm bound one-more Module-ISIS

# Computing the smallest BKZ block size
for beta in range(400,700):
        B = ((beta/(2*math.pi*math.e))*((math.pi*beta)**(1/beta)))**(1/(2*(beta-1)))
        if B**(bkz_d - 1) * q**(1/2) <= bound:
            print("BKZ block size:     ", beta)
            print("root hermite factor:",B)
            print("MSIS security:      ", round(0.292*beta), "bits")
            break

print("Signature size:     ",  round(deg*power/8000,2), "KB")

print("LWE hardness:")
p =  LWE.Parameters(n=k, q=q, Xs=ND.SparseTernary(p=w, m=w), Xe=ND.SparseTernary(p=w/2, m=w/2),m=deg)
LWE.estimate.rough(p)

print("\nFor MSIS scheme:\n")

power = 27
q = 2**power
deg = 1024
k = deg*2        
bkz_d = 2*k 

n_sigma_0 = 1.5*sqrt(q) *2*sqrt(deg)		# norm of sigma_0
s = gamma * n_sigma_0 * (1+sqrt(3*w))		# standard deviation of signature
bound = 4*s*sqrt(deg)				# norm bound for MSIS (2*bound for valid sig)

# Computing the smallest BKZ block size
for beta in range(400,700):
        B = ((beta/(2*math.pi*math.e))*((math.pi*beta)**(1/beta)))**(1/(2*(beta-1)))
        if B**(bkz_d - 1) * q**(1/2) <= bound:
            print("BKZ block size:     ", beta)
            print("root hermite factor:",B)
            print("MSIS security:      ", round(0.292*beta), "bits")
            break

print("Signature size:     ",  round(deg*power/8000,2), "KB")

print("LWE hardness:")
p =  LWE.Parameters(n=k, q=q, Xs=ND.SparseTernary(p=w, m=w), Xe=ND.SparseTernary(p=w/2, m=w/2),m=deg)
LWE.estimate.rough(p)



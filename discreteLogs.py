import os
import sys
import math
from datetime import datetime, timedelta


def main():
    p = input("Input prime for mod value: ")
    l = input("Input log base: ")
    e = input("Input Number: ")
    print("Pick Algorithm:")
    print("[1] Divide and Conquer")
    print("[2] Pohlig-Hellman")
    print("[3] Both")
    algo = input("Enter number: ")
    now = datetime.now()
    start = timedelta(hours=now.hour, minutes=now.minute, seconds=now.second, microseconds=now.microsecond)
    if(algo == 1):
        print("Running Divide and Conquer... \n")
        divideAndConquer(p,l,e)
    elif(algo == 2):
        print("Running Pohlig-Hellman... \n")
        pohlig(p,l,e)
    elif(algo == 3):
        print("Running Pohlig-Hellman... \n")
        pohlig(p,l,e)
        now = datetime.now()
        stop = timedelta(hours=now.hour, minutes=now.minute, seconds=now.second, microseconds=now.microsecond)
        print("\nRan for: "+str(stop-start)+" (H:MM:SS.micros)\n")
        now = datetime.now()
        start = timedelta(hours=now.hour, minutes=now.minute, seconds=now.second, microseconds=now.microsecond)
        print("Running Divide and Conquer... \n")
        divideAndConquer(p,l,e)
    else:
        print("error. not an option")
        raise KeyboardInterrupt
    now = datetime.now()
    stop = timedelta(hours=now.hour, minutes=now.minute, seconds=now.second, microseconds=now.microsecond)
    print("\nRan for: "+str(stop-start)+" (H:MM:SS.micros)")

def pohlig(p,l,e):
    xCongruences = []
    print("Congruences for x:")
    for (fac, pwr) in orderFactors(p-1):
        congruence = getXModP(e,l,p,fac,pwr)
        print("   x mod "+str(congruence[1])+" = "+str(congruence[0]))
        xCongruences.append(congruence)
    x = crt(xCongruences)
    print("For "+str(l)+"^x = "+str(e)+" mod "+str(p)+", x = "+str(x))

def divideAndConquer(p,l,e):
    coef = findEqX(p, e)
    a_values = computeAValues(coef[0],l,p)
    print("a = "+str(coef[0]))
    print(a_values)
    print("")
    b_values = computeBValues(coef[1],l,p, e)
    print("b = "+str(coef[1]))
    print(b_values)
    ab_pair = findPair(a_values, b_values)
    if(ab_pair):
        print("\nFor "+str(l)+"^x = "+str(e)+" mod "+str(p)+":")
        print("a = "+str(coef[0])+", a_i = "+str(ab_pair[0]))
        print("b = "+str(coef[1])+", b_i = "+str(ab_pair[1]))
        print("x = a*(a_i) + b_i")
        print("x = "+str(coef[0])+" * "+str(ab_pair[0])+" + "+str(ab_pair[1])+" = "+str(coef[0]*ab_pair[0] + ab_pair[1]))
    else:
        print("Error, no solution found")

def findEqX(p, e):
    a = ((int)(math.sqrt((p-1))))
    b = a + 1
    return (a,b)

def computeAValues(a, l, p):
    i = 0
    a_values = []
    while(i <= a):
        num = (l ** (a*i))%p
        a_values.append(num)
        i+=1
    return a_values

def computeBValues(b,l,p,e):
    i = 0
    b_values={}
    while(i <= b):
        num = (e*pow(modinv(l,p), i, p))%p
        b_values.update({num: i})
        i+= 1
    return b_values

def findPair(a, b):
    foundPair = ()
    i = 0
    while(i < len(a)):
        try:
            b[a[i]]
            foundPair = (i, b[a[i]])
            break
        except KeyError:
            i += 1
    return foundPair

def egcd(a, b):
    if(a == 0):
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modinv(a, m):
    g, x, y = egcd(a, m)
    if(g != 1):
        raise Exception('modular inverse does not exist')
    else:
        return x % m

def factor(n, startFrom=2):
    if(n <= 1):
        return []
    d = startFrom
    factors = []
    while(n >= d*d):
      if((n % d) == 0):
        factors.append(d)
        n = n/d
      else:
        d += 1 + d % 2
    factors.append(n)
    return factors

def orderFactors(num):
    factors = factor(num)
    if not factors:
        return []
    current = NotImplemented
    n = 0
    pairs = []
    for e in factors:
        if(e == current):
            n += 1
        else:
            if(n > 0):
                pairs.append((current, n))
            n = 1
            current = e
    pairs.append((current, n))
    return pairs

def getXModP(e, l, p, fac, pwr):
    order = (p-1)/fac
    eCurrent = e
    xFinal = 0
    lRaisedModp = pow(l, order, p)
    facPow = 1
    lInv = modinv(l,p)
    for i in range(0,pwr):
        eRaisedModp = pow(eCurrent, order, p)
        xCurrent = bruteForce(lRaisedModp, eRaisedModp, p)
        xFinal += xCurrent*facPow
        eCurrent = eCurrent*pow(lInv, xCurrent*facPow, p) % p
        facPow *= fac
        order /= fac
    return (xFinal,facPow)

def bruteForce(a,b,p):
    a_x = 1
    b %= p
    for x in range(p-1):
        if(a_x == b):
            return x
        a_x = a_x * a % p
    return None

# chinese remainder theorem
def crt(pairs):
    (r1, p1) = pairs[0]
    for (r,p) in pairs[1:]:
        k = ((r-r1)*modinv(p1,p)) % p
        r1 = (r1+p1*k) % (p1*p)
        p1 *= p
    return r1

#Execute the wrapper
if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print 'oops... interrupted \_[o.O]_/'
        try:
            sys.exit(0)
        except SystemExit:
            os._exit(0)

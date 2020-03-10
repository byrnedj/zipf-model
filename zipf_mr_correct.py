import numpy as np
import decimal
import argparse
import math

def compute_normalizer(alpha, m):
    #compute normalizer of Zipf distribution probabilities
    total = decimal.Decimal(0)
    for i in range(1, m+1):
        temp = decimal.Decimal(i)**decimal.Decimal(-alpha)
        total += temp
    norm = decimal.Decimal(1)/total

    #--DEBUG--
    #totalDebug = decimal.Decimal(0)
    #for j in range(1, 25000):
    #    temp = decimal.Decimal(j)**decimal.Decimal(-alpha)
    #    totalDebug += temp
    #    print(j, norm*totalDebug)

    print(norm)
    return norm

def populate_derivatives(m, n):
    #compute derivatives, store in drv
    global drv 
    drv = [decimal.Decimal(0)]*n
    for i in range(0,n,1000):
        for j in range(m):
            drv[i] += (ln_expr[j]*(expr[j]**decimal.Decimal(i)))

def populate_subexpressions(m, alpha, norm):
    #compute reused subexpressions (2^.. terms). **optimization**
    global expr 
    global ln_expr
    ln_expr = [decimal.Decimal(0)]*m
    expr = [decimal.Decimal(1) - (decimal.Decimal(norm)*(decimal.Decimal(i)**decimal.Decimal(-alpha))) for i in range(1, m+1)]
    #for x in range(len(expr)):
        #print(x, expr[x])
    #    ln_expr.append(decimal.Decimal(-1)*expr[x].ln())
    ln_expr = [decimal.Decimal(-1)*x.ln() for x in expr]

def populate_footprint(m, n):
    #compute footprints, store in fp
    global fp 
    fp = [decimal.Decimal(0)]*n
    for i in range(0,n,1000):
        for j in range(m):
            fp[i] += (decimal.Decimal(1) - (expr[j]**decimal.Decimal(i)))


def compute_mrs(args):
    #populate arrays, compute normalizer
    norm = compute_normalizer(float(args.alpha), int(args.m))
    populate_subexpressions(int(args.m), float(args.alpha), norm)
    populate_footprint(int(args.m), int(args.n))
    populate_derivatives(int(args.m), int(args.n))
    
    #print parameter info
    print("parameters: data size =", args.m, "trace length =", args.n, "zipf parameter =", args.alpha)

    #print miss ratios predicted by both continuous and discrete analyses of footprint values
    cache_size = 0
    for x in range(int(args.n)-1):
        #only print MR one time per whole-valued cache size
        if math.floor(fp[x]) > cache_size:
            print('cache size:', fp[x], '\n    miss ratio (discrete):', fp[x+1] - fp[x], '\n    miss ratio (continuous):', drv[x])
            cache_size = math.floor(fp[x])


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Zipfian distribution miss ratio generator")
    p.add_argument("m", help="data size")
    p.add_argument("n", help="trace length")
    p.add_argument("alpha", help="zipf distribution parameter")
    decimal.getcontext().prec = 100
    args = p.parse_args()
    compute_mrs(args)
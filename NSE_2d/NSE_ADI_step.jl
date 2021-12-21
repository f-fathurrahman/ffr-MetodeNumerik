# modifies fi
function NSE_ADI_step(ami, api, alph, xs2, fi_)

    fi = copy(fi_)' # transpose
    m = size(fi,1)
    n = size(fi,2)

    # Forward substitution
    for i in 2:n
        @views fi[:,i] = fi[:,i] - ami[i]*fi[:,i-1]
    end

    # Backward substitution
    fi[:,n] = fi[:,n]*alph[n]
    for i in n-1:-1:1
        @views fi[:,i] = fi[:,i]*alph[i] - api[i]*fi[:,i+1]
    end

    # Final solution
    #xs1 = zeros(m,1);
    xs1 = ( fi[:,1] + fi[:,n] ) / (1.0 + xs2[1] + xs2[n])
    for i in 1:n
        @views fi[:,i] = fi[:,i] - xs1*xs2[i]
    end

    return fi
end
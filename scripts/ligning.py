hus_pris = 2.8*10**6
egenkapital = 0.2*hus_pris
debitorrente = 0.025
fradrag = 0.336
maks_rente_udgift = 100000

lån = hus_pris-egenkapital

print (egenkapital)
print (lån)

rente_udgift = lån*debitorrente

print (rente_udgift)

print (rente_udgift*1.25 + 3200*12)
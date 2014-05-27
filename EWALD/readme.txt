zdrojak s r- i k-space Ewaldem a testovacim main():
  ewaldf.c

include moduly:
  ertest.c : testovaci tisk splajnu
  ewaldf.h : hlavicky, struktury
  ewaldprt.c : tisky k-space
  ewpass1.c, ewpass2.c , ewpass2m.c : k-space casti
  include.h : MACSIMUS interface 
  loop.h : smycky
  mystring.c : retezce
  vector3d.h : makra pro vektory

exe:
  ewaldf : linux 64 bit exe
  ewald : primitivni testovaci program, linux 64 bit exe

testovaci data:
  ewald.dat : 5 naboju
  1000.dat : Na500Cl500 melt


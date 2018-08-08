SHELL = /bin/sh

OBJS_f90d1 = OXYG6.o BRDFGRID.o SAND.o METH2.o MOCA5.o DMC.o OZON1.o BBM.o SCATRA.o \
METH3.o MOCA6.o DUST.o PLANPOL.o AEROPROF.o SEAWIFS.o METH4.o MODIS.o ENVIRO.o POLDER.o \
CAVIS.o SOLIRR.o METH5.o MODIS1km.o EQUIVWL.o POLGLIT.o CHAND.o SOOT.o METH6.o MODISALBE.o \
ETM.o POLNAD.o ARTALBE.o SPECINTERP.o MIDSUM.o MODISBRDF.o GAUSS.o POSGE.o ALI.o \
SPLIE2.o MIDWIN.o MODISBRDFFOS.o GLI.o POSGW.o ATMREF.o SPLIN2.o MIE.o MSS.o GOES.o POSLAN.o \
ABSTRA.o SPLINE.o NIOX1.o HAPKALBE.o POSMTO.o AATSR.o SPLINT.o NIOX2.o HAPKBRDF.o POSNOA.o \
STM.o NIOX3.o VERSBRDF.o HRV.o POSSOL.o SUBSUM.o NIOX4.o VERSTOOLS.o HYPBLUE.o POSSPO.o SUBWIN.o \
NIOX5.o VGT.o IAPIALBE.o PRESPLANE.o TM.o NIOX6.o VIIRS.o PRESSURE.o TROPIC.o OCEA.o WALTALBE.o \
PRINT_ERROR.o TRUNCA.o OCEAALBE.o AVHRR.o WALTBRDF.o PROSAILALBE.o IAPIBRDF.o US62.o OCEABRDF.o \
ASTER.o WATE.o PROSAILBRDF.o IAPITOOLS.o VARSOL.o OCEABRDFFAST.o ACRMBRDF.o WAVA1.o PROSAILBRDFFOS.o \
INTERP.o VEGETA.o OCEATOOLS.o BDM.o WAVA2.o PROSAILTOOLS.o ISO.o VERSALBE.o ODA550.o ACRMBRDFFOS.o \
WAVA3.o RAHMALBE.o KERNEL.o ODRAYL.o ACRMALBE.o WAVA4.o RAHMBRDF.o KERNELPOL.o CLEARW.o OS.o \
ARTBRDFFOS.o WAVA5.o RHAMBRDFFOS.o LAKEW.o MINNALBE.o CSALBR.o OSPOL.o ARTTOOL.o WAVA6.o \
ROSSLIMAIGNANALBE.o LDCM.o MINNBRDF.o DICA1.o OSSURF.o AEROSO.o ROSSLIMAIGNANBRDF.o MAS.o \
MOCA1.o DICA2.o OXYG3.o ACRMTOOL.o paramdef.o ROSSLIMAIGNANBRDFFOS.o MERIS.o MOCA2.o DICA3.o OXYG4.o \
AKTOOL.o ROUJALBE.o METEO.o MOCA3.o DISCOM.o OXYG5.o ARTBRDF.o ROUJBRDF.o METH1.o MOCA4.o DISCRE.o

OBJS_f90d2 = main.o
OBJS_f90d3 = mainlutaero.o

EXE2 = sixs90
EXE3 = sixslutaero90
FC = gfortran

CFLAGS = -Wall -g -Wtabs -ffpe-trap=zero,overflow,invalid,denormal,underflow 

all : $(EXE2) $(EXE3)

$(EXE2) : $(OBJS_f90d1) $(OBJS_f90d2)
	$(FC) -o $(EXE2) $(OBJS_f90d2) $(OBJS_f90d1)

$(EXE3) : $(OBJS_f90d1) $(OBJS_f90d3)
	$(FC) -o $(EXE3) $(OBJS_f90d3) $(OBJS_f90d1)

$(OBJS_f90d1) $(OBJS_f90d2) $(OBJS_f90d3):
	$(FC) $(CFLAGS) -c $(@:.o=.f90) -o $@


clean :
	rm -f *.o
	rm -f *.mod

# Dependencies of files
mainlutaero.o: \
    mainlutaero.f90 \
    paramdef.o
OXYG6.o: \
    OXYG6.f90
BRDFGRID.o: \
    BRDFGRID.f90
SAND.o: \
    SAND.f90
METH2.o: \
    METH2.f90
MOCA5.o: \
    MOCA5.f90
DMC.o: \
    DMC.f90
OZON1.o: \
    OZON1.f90
BBM.o: \
    BBM.f90 \
    paramdef.o
SCATRA.o: \
    SCATRA.f90
METH3.o: \
    METH3.f90
MOCA6.o: \
    MOCA6.f90
DUST.o: \
    DUST.f90 \
    paramdef.o
PLANPOL.o: \
    PLANPOL.f90
AEROPROF.o: \
    AEROPROF.f90 \
    paramdef.o
SEAWIFS.o: \
    SEAWIFS.f90
METH4.o: \
    METH4.f90
MODIS.o: \
    MODIS.f90
ENVIRO.o: \
    ENVIRO.f90
POLDER.o: \
    POLDER.f90
CAVIS.o: \
    CAVIS.f90
SOLIRR.o: \
    SOLIRR.f90
METH5.o: \
    METH5.f90
MODIS1km.o: \
    MODIS1km.f90
EQUIVWL.o: \
    EQUIVWL.f90
POLGLIT.o: \
    POLGLIT.f90
CHAND.o: \
    CHAND.f90
SOOT.o: \
    SOOT.f90 \
    paramdef.o
METH6.o: \
    METH6.f90
MODISALBE.o: \
    MODISALBE.f90
ETM.o: \
    ETM.f90
POLNAD.o: \
    POLNAD.f90
ARTALBE.o: \
    ARTALBE.f90
SPECINTERP.o: \
    SPECINTERP.f90 \
    paramdef.o
MIDSUM.o: \
    MIDSUM.f90
MODISBRDF.o: \
    MODISBRDF.f90
GAUSS.o: \
    GAUSS.f90
POSGE.o: \
    POSGE.f90
ALI.o: \
    ALI.f90
SPLIE2.o: \
    SPLIE2.f90
MIDWIN.o: \
    MIDWIN.f90
MODISBRDFFOS.o: \
    MODISBRDFFOS.f90 \
    paramdef.o
GLI.o: \
    GLI.f90
POSGW.o: \
    POSGW.f90
ATMREF.o: \
    ATMREF.f90
SPLIN2.o: \
    SPLIN2.f90
MIE.o: \
    MIE.f90 \
    paramdef.o
MSS.o: \
    MSS.f90
GOES.o: \
    GOES.f90
POSLAN.o: \
    POSLAN.f90
ABSTRA.o: \
    ABSTRA.f90
SPLINE.o: \
    SPLINE.f90
NIOX1.o: \
    NIOX1.f90
HAPKALBE.o: \
    HAPKALBE.f90
POSMTO.o: \
    POSMTO.f90
AATSR.o: \
    AATSR.f90
SPLINT.o: \
    SPLINT.f90
NIOX2.o: \
    NIOX2.f90
HAPKBRDF.o: \
    HAPKBRDF.f90
POSNOA.o: \
    POSNOA.f90
STM.o: \
    STM.f90 \
    paramdef.o
NIOX3.o: \
    NIOX3.f90
VERSBRDF.o: \
    VERSBRDF.f90
HRV.o: \
    HRV.f90
POSSOL.o: \
    POSSOL.f90
SUBSUM.o: \
    SUBSUM.f90
NIOX4.o: \
    NIOX4.f90
VERSTOOLS.o: \
    VERSTOOLS.f90
HYPBLUE.o: \
    HYPBLUE.f90
POSSPO.o: \
    POSSPO.f90
SUBWIN.o: \
    SUBWIN.f90
NIOX5.o: \
    NIOX5.f90
VGT.o: \
    VGT.f90
IAPIALBE.o: \
    IAPIALBE.f90
PRESPLANE.o: \
    PRESPLANE.f90
TM.o: \
    TM.f90
NIOX6.o: \
    NIOX6.f90
VIIRS.o: \
    VIIRS.f90
PRESSURE.o: \
    PRESSURE.f90
TROPIC.o: \
    TROPIC.f90
OCEA.o: \
    OCEA.f90 \
    paramdef.o
WALTALBE.o: \
    WALTALBE.f90
PRINT_ERROR.o: \
    PRINT_ERROR.f90
TRUNCA.o: \
    TRUNCA.f90 \
    paramdef.o
OCEAALBE.o: \
    OCEAALBE.f90
AVHRR.o: \
    AVHRR.f90
WALTBRDF.o: \
    WALTBRDF.f90
PROSAILALBE.o: \
    PROSAILALBE.f90
IAPIBRDF.o: \
    IAPIBRDF.f90
US62.o: \
    US62.f90
OCEABRDF.o: \
    OCEABRDF.f90
ASTER.o: \
    ASTER.f90
WATE.o: \
    WATE.f90 \
    paramdef.o
PROSAILBRDF.o: \
    PROSAILBRDF.f90
IAPITOOLS.o: \
    IAPITOOLS.f90
VARSOL.o: \
    VARSOL.f90
OCEABRDFFAST.o: \
    OCEABRDFFAST.f90
ACRMBRDF.o: \
    ACRMBRDF.f90
WAVA1.o: \
    WAVA1.f90
PROSAILBRDFFOS.o: \
    PROSAILBRDFFOS.f90 \
    paramdef.o
INTERP.o: \
    INTERP.f90 \
    paramdef.o
VEGETA.o: \
    VEGETA.f90
OCEATOOLS.o: \
    OCEATOOLS.f90
BDM.o: \
    BDM.f90 \
    paramdef.o
WAVA2.o: \
    WAVA2.f90
PROSAILTOOLS.o: \
    PROSAILTOOLS.f90
ISO.o: \
    ISO.f90 \
    paramdef.o
VERSALBE.o: \
    VERSALBE.f90
ODA550.o: \
    ODA550.f90
ACRMBRDFFOS.o: \
    ACRMBRDFFOS.f90 \
    paramdef.o
WAVA3.o: \
    WAVA3.f90
RAHMALBE.o: \
    RAHMALBE.f90
KERNEL.o: \
    KERNEL.f90 \
    paramdef.o
ODRAYL.o: \
    ODRAYL.f90
ACRMALBE.o: \
    ACRMALBE.f90
WAVA4.o: \
    WAVA4.f90
RAHMBRDF.o: \
    RAHMBRDF.f90
KERNELPOL.o: \
    KERNELPOL.f90 \
    paramdef.o
CLEARW.o: \
    CLEARW.f90
OS.o: \
    OS.f90 \
    paramdef.o
ARTBRDFFOS.o: \
    ARTBRDFFOS.f90 \
    paramdef.o
WAVA5.o: \
    WAVA5.f90
RHAMBRDFFOS.o: \
    RHAMBRDFFOS.f90 \
    paramdef.o
LAKEW.o: \
    LAKEW.f90
MINNALBE.o: \
    MINNALBE.f90
CSALBR.o: \
    CSALBR.f90
OSPOL.o: \
    OSPOL.f90 \
    paramdef.o
ARTTOOL.o: \
    ARTTOOL.f90
WAVA6.o: \
    WAVA6.f90
ROSSLIMAIGNANALBE.o: \
    ROSSLIMAIGNANALBE.f90
LDCM.o: \
    LDCM.f90
MINNBRDF.o: \
    MINNBRDF.f90
DICA1.o: \
    DICA1.f90
OSSURF.o: \
    OSSURF.f90 \
    paramdef.o
AEROSO.o: \
    AEROSO.f90 \
    paramdef.o
main.o: \
    main.f90 \
    paramdef.o
ROSSLIMAIGNANBRDF.o: \
    ROSSLIMAIGNANBRDF.f90
MAS.o: \
    MAS.f90
MOCA1.o: \
    MOCA1.f90
DICA2.o: \
    DICA2.f90
OXYG3.o: \
    OXYG3.f90
ACRMTOOL.o: \
    ACRMTOOL.f90
paramdef.o: \
    paramdef.f90
ROSSLIMAIGNANBRDFFOS.o: \
    ROSSLIMAIGNANBRDFFOS.f90 \
    paramdef.o
MERIS.o: \
    MERIS.f90
MOCA2.o: \
    MOCA2.f90
DICA3.o: \
    DICA3.f90
OXYG4.o: \
    OXYG4.f90
AKTOOL.o: \
    AKTOOL.f90
ROUJALBE.o: \
    ROUJALBE.f90
METEO.o: \
    METEO.f90
MOCA3.o: \
    MOCA3.f90
DISCOM.o: \
    DISCOM.f90 \
    paramdef.o
OXYG5.o: \
    OXYG5.f90
ARTBRDF.o: \
    ARTBRDF.f90
ROUJBRDF.o: \
    ROUJBRDF.f90
METH1.o: \
    METH1.f90
MOCA4.o: \
    MOCA4.f90
DISCRE.o: \
    DISCRE.f90

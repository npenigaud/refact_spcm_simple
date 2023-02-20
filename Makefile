include Makefile.inc

./abor1.intfb.ok: ./abor1.intfb.h 
	touch ./abor1.intfb.ok

./mxmaop.intfb.ok: ./mxmaop.intfb.h 
	touch ./mxmaop.intfb.ok

./mxptma.intfb.ok: ./mxptma.intfb.h 
	touch ./mxptma.intfb.ok

./mxture.intfb.ok: ./mxture.intfb.h 
	touch ./mxture.intfb.ok

./mxturs.intfb.ok: ./mxturs.intfb.h 
	touch ./mxturs.intfb.ok

./set2pe.intfb.ok: ./set2pe.intfb.h yommp0.o yomlun.o
	touch ./set2pe.intfb.ok

./sgemmx.intfb.ok: ./sgemmx.intfb.h 
	touch ./sgemmx.intfb.ok

./sigam_gp.intfb.ok: ./sigam_gp.intfb.h geometry_mod.o yomcst.o yomdyn.o
	touch ./sigam_gp.intfb.ok

./sigam_sp_openmp.intfb.ok: ./sigam_sp_openmp.intfb.h geometry_mod.o yomcst.o yomdyn.o
	touch ./sigam_sp_openmp.intfb.ok

./simplico.intfb.ok: ./simplico.intfb.h 
	touch ./simplico.intfb.ok

./sitnu_gp.intfb.ok: ./sitnu_gp.intfb.h geometry_mod.o yomcst.o yomdyn.o
	touch ./sitnu_gp.intfb.ok

./sitnu_sp_openmp.intfb.ok: ./sitnu_sp_openmp.intfb.h geometry_mod.o yomcst.o yomdyn.o
	touch ./sitnu_sp_openmp.intfb.ok

./spcimpfsolve.intfb.ok: ./spcimpfsolve.intfb.h geometry_mod.o yomdyn.o yomrip.o yomcst.o
	touch ./spcimpfsolve.intfb.ok

./spcm_simple.intfb.ok: ./spcm_simple.intfb.h type_model.o geometry_mod.o
	touch ./spcm_simple.intfb.ok

./spcsi.intfb.ok: ./spcsi.intfb.h geometry_mod.o yommp0.o yomdyn.o yomlddh.o yomrip.o yomcst.o
	touch ./spcsi.intfb.ok

./trmtos.intfb.ok: ./trmtos.intfb.h geometry_mod.o yommp0.o yomtag.o
	touch ./trmtos.intfb.ok

./trstom.intfb.ok: ./trstom.intfb.h geometry_mod.o yommp0.o yomtag.o
	touch ./trstom.intfb.ok

./verder.intfb.ok: ./verder.intfb.h 
	touch ./verder.intfb.ok

./verdisint.intfb.ok: ./verdisint.intfb.h yomcver.o yomlun.o yomvert.o
	touch ./verdisint.intfb.ok

./verint.intfb.ok: ./verint.intfb.h yomlun.o
	touch ./verint.intfb.ok

mxmaop.ok: mxmaop.h
	touch mxmaop.ok
mxptma.ok: mxptma.h
	touch mxptma.ok
mxture.ok: mxture.h
	touch mxture.ok
mxturs.ok: mxturs.h
	touch mxturs.ok
simplico.ok: simplico.h
	touch simplico.ok
abor1.o: abor1.F90 
	$(FC) -c abor1.F90

geometry_mod.o: geometry_mod.F90 type_geometry.o yomlun.o
	$(FC) -c geometry_mod.F90

intdyn_mod.o: intdyn_mod.F90 
	$(FC) -c intdyn_mod.F90

model_diagnostics_mod.o: model_diagnostics_mod.F90 yomlddh.o yomspddh.o
	$(FC) -c model_diagnostics_mod.F90

model_dynamics_mod.o: model_dynamics_mod.F90 yomdyna.o yomdyn.o
	$(FC) -c model_dynamics_mod.F90

model_general_conf_mod.o: model_general_conf_mod.F90 yomrip.o
	$(FC) -c model_general_conf_mod.F90

mxmaop.o: mxmaop.F90 
	$(FC) -c mxmaop.F90

mxptma.o: mxptma.F90 
	$(FC) -c mxptma.F90

mxture.o: mxture.F90 
	$(FC) -c mxture.F90

mxturs.o: mxturs.F90 mxture.ok
	$(FC) -c mxturs.F90

myrecvset.o: myrecvset.F90 abor1.intfb.ok
	$(FC) -c myrecvset.F90

mysendset.o: mysendset.F90 abor1.intfb.ok
	$(FC) -c mysendset.F90

reglatlon_field_mix.o: reglatlon_field_mix.F90 yomlun.o
	$(FC) -c reglatlon_field_mix.F90

set2pe.o: set2pe.F90 yommp0.o yomlun.o abor1.intfb.ok
	$(FC) -c set2pe.F90

sgemmx.o: sgemmx.F90 
	$(FC) -c sgemmx.F90

sigam_gp.o: sigam_gp.F90 geometry_mod.o yomcst.o yomdyn.o verdisint.intfb.ok
	$(FC) -c sigam_gp.F90

sigam_sp_openmp.o: sigam_sp_openmp.F90 geometry_mod.o yomcst.o yomdyn.o verdisint.intfb.ok
	$(FC) -c sigam_sp_openmp.F90

simplico.o: simplico.F90 
	$(FC) -c simplico.F90

sitnu_gp.o: sitnu_gp.F90 geometry_mod.o yomcst.o yomdyn.o verdisint.intfb.ok
	$(FC) -c sitnu_gp.F90

sitnu_sp_openmp.o: sitnu_sp_openmp.F90 geometry_mod.o yomcst.o yomdyn.o verdisint.intfb.ok
	$(FC) -c sitnu_sp_openmp.F90

spcimpfsolve.o: spcimpfsolve.F90 geometry_mod.o yomdyn.o yomrip.o yomcst.o simplico.ok trmtos.intfb.ok trstom.intfb.ok
	$(FC) -c spcimpfsolve.F90

spcm.o: spcm.F90 xrd_getoptions.o type_model.o geometry_mod.o util_model_mod.o util_geometry_mod.o util_yommp0_mod.o spcm_simple.intfb.ok abor1.intfb.ok
	$(FC) -c spcm.F90

spcm_simple.o: spcm_simple.F90 type_model.o geometry_mod.o spcsi.intfb.ok trmtos.intfb.ok trstom.intfb.ok
	$(FC) -c spcm_simple.F90

spcsi.o: spcsi.F90 geometry_mod.o yommp0.o yomdyn.o yomlddh.o yomrip.o yomcst.o mxmaop.ok mxptma.ok mxture.ok mxturs.ok abor1.intfb.ok sigam_sp_openmp.intfb.ok spcimpfsolve.intfb.ok sitnu_sp_openmp.intfb.ok
	$(FC) -c spcsi.F90

trmtos.o: trmtos.F90 geometry_mod.o yommp0.o yomtag.o abor1.intfb.ok set2pe.intfb.ok
	$(FC) -c trmtos.F90

trstom.o: trstom.F90 geometry_mod.o yommp0.o yomtag.o abor1.intfb.ok set2pe.intfb.ok
	$(FC) -c trstom.F90

type_geometry.o: type_geometry.F90 yomvert.o yomsta.o yomlap.o yomleg.o yomdim.o yomdimv.o yommp.o yomgem.o yomcsgeom.o yomgsgeom.o yomorog.o type_spgeom.o yemdim.o yemgeo.o yemmp.o yemlap.o yemgsl.o yemlbc_geo.o yomcver.o
	$(FC) -c type_geometry.F90

type_model.o: type_model.F90 model_general_conf_mod.o model_dynamics_mod.o model_diagnostics_mod.o yomcst.o
	$(FC) -c type_model.F90

type_spgeom.o: type_spgeom.F90 
	$(FC) -c type_spgeom.F90

util_geometry_mod.o: util_geometry_mod.F90 type_geometry.o util_tcsgeom_blocked_mod.o util_tcsgeom_mod.o util_tcsgleg_mod.o util_tcver_mod.o util_tdimv_mod.o util_tdim_mod.o util_tedim_mod.o util_tegeo_mod.o util_tegsl_mod.o util_telbc_geo_mod.o util_temmp_mod.o util_tgem_mod.o util_tgsgeom_blocked_mod.o util_tgsgeom_mod.o util_tlap_mod.o util_tlep_mod.o util_tmp_mod.o util_torog_blocked_mod.o util_torog_mod.o util_tspgeom_mod.o util_tsta_mod.o util_tvab_mod.o util_tvertical_geom_mod.o util_tveta_mod.o util_tvfe_mod.o
	$(FC) -c util_geometry_mod.F90

util_model_diagnostics_type_mod.o: util_model_diagnostics_type_mod.F90 model_diagnostics_mod.o util_tlddh_mod.o util_tspddh_mod.o
	$(FC) -c util_model_diagnostics_type_mod.F90

util_model_dynamics_type_mod.o: util_model_dynamics_type_mod.F90 model_dynamics_mod.o util_tdyna_mod.o util_tdyn_mod.o
	$(FC) -c util_model_dynamics_type_mod.F90

util_model_general_conf_type_mod.o: util_model_general_conf_type_mod.F90 model_general_conf_mod.o util_trip_mod.o
	$(FC) -c util_model_general_conf_type_mod.F90

util_model_mod.o: util_model_mod.F90 type_model.o util_model_diagnostics_type_mod.o util_model_dynamics_type_mod.o util_model_general_conf_type_mod.o util_tcst_mod.o
	$(FC) -c util_model_mod.F90

util_reglatlon_field_mod.o: util_reglatlon_field_mod.F90 reglatlon_field_mix.o
	$(FC) -c util_reglatlon_field_mod.F90

util_tcsgeom_blocked_mod.o: util_tcsgeom_blocked_mod.F90 yomcsgeom.o
	$(FC) -c util_tcsgeom_blocked_mod.F90

util_tcsgeom_mod.o: util_tcsgeom_mod.F90 yomcsgeom.o
	$(FC) -c util_tcsgeom_mod.F90

util_tcsgleg_mod.o: util_tcsgleg_mod.F90 yomleg.o
	$(FC) -c util_tcsgleg_mod.F90

util_tcst_mod.o: util_tcst_mod.F90 yomcst.o
	$(FC) -c util_tcst_mod.F90

util_tcty_mod.o: util_tcty_mod.F90 intdyn_mod.o
	$(FC) -c util_tcty_mod.F90

util_tcver_mod.o: util_tcver_mod.F90 yomcver.o
	$(FC) -c util_tcver_mod.F90

util_tdim_mod.o: util_tdim_mod.F90 yomdim.o
	$(FC) -c util_tdim_mod.F90

util_tdimv_mod.o: util_tdimv_mod.F90 yomdimv.o
	$(FC) -c util_tdimv_mod.F90

util_tdyn_mod.o: util_tdyn_mod.F90 yomdyn.o
	$(FC) -c util_tdyn_mod.F90

util_tdyna_mod.o: util_tdyna_mod.F90 yomdyna.o util_tgflt_mod.o util_tgmvt_mod.o util_ttnd_mod.o
	$(FC) -c util_tdyna_mod.F90

util_teaerc_macc_mod.o: util_teaerc_macc_mod.F90 yoeaerc.o
	$(FC) -c util_teaerc_macc_mod.F90

util_teaerc_tegen_mod.o: util_teaerc_tegen_mod.F90 yoeaerc.o
	$(FC) -c util_teaerc_tegen_mod.F90

util_tecmip_mod.o: util_tecmip_mod.F90 yoecmip.o
	$(FC) -c util_tecmip_mod.F90

util_tedim_mod.o: util_tedim_mod.F90 yemdim.o
	$(FC) -c util_tedim_mod.F90

util_tegeo_mod.o: util_tegeo_mod.F90 yemgeo.o
	$(FC) -c util_tegeo_mod.F90

util_tegsl_mod.o: util_tegsl_mod.F90 yemgsl.o
	$(FC) -c util_tegsl_mod.F90

util_telbc_geo_mod.o: util_telbc_geo_mod.F90 yemlbc_geo.o
	$(FC) -c util_telbc_geo_mod.F90

util_temmp_mod.o: util_temmp_mod.F90 yemmp.o
	$(FC) -c util_temmp_mod.F90

util_teozoc_mod.o: util_teozoc_mod.F90 yoeozoc.o
	$(FC) -c util_teozoc_mod.F90

util_tgem_mod.o: util_tgem_mod.F90 yomgem.o
	$(FC) -c util_tgem_mod.F90

util_tgflt_mod.o: util_tgflt_mod.F90 intdyn_mod.o
	$(FC) -c util_tgflt_mod.F90

util_tgmvt_mod.o: util_tgmvt_mod.F90 intdyn_mod.o
	$(FC) -c util_tgmvt_mod.F90

util_tgsgeom_blocked_mod.o: util_tgsgeom_blocked_mod.F90 yomgsgeom.o
	$(FC) -c util_tgsgeom_blocked_mod.F90

util_tgsgeom_mod.o: util_tgsgeom_mod.F90 yomgsgeom.o
	$(FC) -c util_tgsgeom_mod.F90

util_thwind_mod.o: util_thwind_mod.F90 intdyn_mod.o
	$(FC) -c util_thwind_mod.F90

util_tlap_mod.o: util_tlap_mod.F90 yomlap.o
	$(FC) -c util_tlap_mod.F90

util_tlddh_mod.o: util_tlddh_mod.F90 yomlddh.o
	$(FC) -c util_tlddh_mod.F90

util_tlep_mod.o: util_tlep_mod.F90 yemlap.o
	$(FC) -c util_tlep_mod.F90

util_tmp_mod.o: util_tmp_mod.F90 yommp.o
	$(FC) -c util_tmp_mod.F90

util_torog_blocked_mod.o: util_torog_blocked_mod.F90 yomorog.o
	$(FC) -c util_torog_blocked_mod.F90

util_torog_mod.o: util_torog_mod.F90 yomorog.o
	$(FC) -c util_torog_mod.F90

util_tpg_type_mod.o: util_tpg_type_mod.F90 intdyn_mod.o
	$(FC) -c util_tpg_type_mod.F90

util_tradghg_mod.o: util_tradghg_mod.F90 yoeradghg.o
	$(FC) -c util_tradghg_mod.F90

util_trcp_mod.o: util_trcp_mod.F90 intdyn_mod.o
	$(FC) -c util_trcp_mod.F90

util_trip_mod.o: util_trip_mod.F90 yomrip.o util_reglatlon_field_mod.o util_teaerc_macc_mod.o util_teaerc_tegen_mod.o util_tecmip_mod.o util_teozoc_mod.o util_tradghg_mod.o
	$(FC) -c util_trip_mod.F90

util_tspddh_mod.o: util_tspddh_mod.F90 yomspddh.o
	$(FC) -c util_tspddh_mod.F90

util_tspgeom_mod.o: util_tspgeom_mod.F90 type_spgeom.o
	$(FC) -c util_tspgeom_mod.F90

util_tsta_mod.o: util_tsta_mod.F90 yomsta.o
	$(FC) -c util_tsta_mod.F90

util_ttnd_mod.o: util_ttnd_mod.F90 intdyn_mod.o
	$(FC) -c util_ttnd_mod.F90

util_tvab_mod.o: util_tvab_mod.F90 yomvert.o
	$(FC) -c util_tvab_mod.F90

util_tvertical_geom_mod.o: util_tvertical_geom_mod.F90 yomvert.o util_tcver_mod.o util_tvab_mod.o util_tveta_mod.o util_tvfe_mod.o
	$(FC) -c util_tvertical_geom_mod.F90

util_tveta_mod.o: util_tveta_mod.F90 yomvert.o
	$(FC) -c util_tveta_mod.F90

util_tvfe_mod.o: util_tvfe_mod.F90 yomvert.o
	$(FC) -c util_tvfe_mod.F90

util_txyb_mod.o: util_txyb_mod.F90 intdyn_mod.o
	$(FC) -c util_txyb_mod.F90

util_txybder_mod.o: util_txybder_mod.F90 intdyn_mod.o
	$(FC) -c util_txybder_mod.F90

util_yommp0_mod.o: util_yommp0_mod.F90 yommp0.o
	$(FC) -c util_yommp0_mod.F90

verder.o: verder.F90 
	$(FC) -c verder.F90

verdisint.o: verdisint.F90 yomcver.o yomlun.o yomvert.o abor1.intfb.ok verder.intfb.ok verint.intfb.ok
	$(FC) -c verdisint.F90

verint.o: verint.F90 yomlun.o abor1.intfb.ok
	$(FC) -c verint.F90

xrd_getoptions.o: xrd_getoptions.F90 xrd_unix_env.o
	$(FC) -c xrd_getoptions.F90

xrd_unix_env.o: xrd_unix_env.F90 
	$(FC) -c xrd_unix_env.F90

yemdim.o: yemdim.F90 
	$(FC) -c yemdim.F90

yemgeo.o: yemgeo.F90 
	$(FC) -c yemgeo.F90

yemgsl.o: yemgsl.F90 
	$(FC) -c yemgsl.F90

yemlap.o: yemlap.F90 
	$(FC) -c yemlap.F90

yemlbc_geo.o: yemlbc_geo.F90 
	$(FC) -c yemlbc_geo.F90

yemmp.o: yemmp.F90 
	$(FC) -c yemmp.F90

yoeaerc.o: yoeaerc.F90 
	$(FC) -c yoeaerc.F90

yoecmip.o: yoecmip.F90 
	$(FC) -c yoecmip.F90

yoeozoc.o: yoeozoc.F90 
	$(FC) -c yoeozoc.F90

yoeradghg.o: yoeradghg.F90 yomcst.o
	$(FC) -c yoeradghg.F90

yomcsgeom.o: yomcsgeom.F90 
	$(FC) -c yomcsgeom.F90

yomcst.o: yomcst.F90 
	$(FC) -c yomcst.F90

yomcver.o: yomcver.F90 yomlun.o
	$(FC) -c yomcver.F90

yomdim.o: yomdim.F90 
	$(FC) -c yomdim.F90

yomdimv.o: yomdimv.F90 
	$(FC) -c yomdimv.F90

yomdyn.o: yomdyn.F90 
	$(FC) -c yomdyn.F90

yomdyna.o: yomdyna.F90 intdyn_mod.o
	$(FC) -c yomdyna.F90

yomgem.o: yomgem.F90 
	$(FC) -c yomgem.F90

yomgsgeom.o: yomgsgeom.F90 
	$(FC) -c yomgsgeom.F90

yomlap.o: yomlap.F90 
	$(FC) -c yomlap.F90

yomlddh.o: yomlddh.F90 
	$(FC) -c yomlddh.F90

yomleg.o: yomleg.F90 
	$(FC) -c yomleg.F90

yomlun.o: yomlun.F90 yomlun_ifsaux.o
	$(FC) -c yomlun.F90

yomlun_ifsaux.o: yomlun_ifsaux.F90 
	$(FC) -c yomlun_ifsaux.F90

yommp.o: yommp.F90 
	$(FC) -c yommp.F90

yommp0.o: yommp0.F90 
	$(FC) -c yommp0.F90

yomorog.o: yomorog.F90 
	$(FC) -c yomorog.F90

yomrip.o: yomrip.F90 yoeozoc.o yoecmip.o yoeradghg.o yoeaerc.o reglatlon_field_mix.o
	$(FC) -c yomrip.F90

yomspddh.o: yomspddh.F90 
	$(FC) -c yomspddh.F90

yomsta.o: yomsta.F90 
	$(FC) -c yomsta.F90

yomtag.o: yomtag.F90 
	$(FC) -c yomtag.F90

yomvert.o: yomvert.F90 yomcver.o
	$(FC) -c yomvert.F90

spcm.x: spcm.o abor1.o geometry_mod.o intdyn_mod.o model_diagnostics_mod.o model_dynamics_mod.o model_general_conf_mod.o mxmaop.o mxptma.o mxture.o mxturs.o myrecvset.o mysendset.o reglatlon_field_mix.o set2pe.o sgemmx.o sigam_gp.o sigam_sp_openmp.o simplico.o sitnu_gp.o sitnu_sp_openmp.o spcimpfsolve.o spcm_simple.o spcsi.o trmtos.o trstom.o type_geometry.o type_model.o type_spgeom.o util_geometry_mod.o util_model_diagnostics_type_mod.o util_model_dynamics_type_mod.o util_model_general_conf_type_mod.o util_model_mod.o util_reglatlon_field_mod.o util_tcsgeom_blocked_mod.o util_tcsgeom_mod.o util_tcsgleg_mod.o util_tcst_mod.o util_tcty_mod.o util_tcver_mod.o util_tdim_mod.o util_tdimv_mod.o util_tdyn_mod.o util_tdyna_mod.o util_teaerc_macc_mod.o util_teaerc_tegen_mod.o util_tecmip_mod.o util_tedim_mod.o util_tegeo_mod.o util_tegsl_mod.o util_telbc_geo_mod.o util_temmp_mod.o util_teozoc_mod.o util_tgem_mod.o util_tgflt_mod.o util_tgmvt_mod.o util_tgsgeom_blocked_mod.o util_tgsgeom_mod.o util_thwind_mod.o util_tlap_mod.o util_tlddh_mod.o util_tlep_mod.o util_tmp_mod.o util_torog_blocked_mod.o util_torog_mod.o util_tpg_type_mod.o util_tradghg_mod.o util_trcp_mod.o util_trip_mod.o util_tspddh_mod.o util_tspgeom_mod.o util_tsta_mod.o util_ttnd_mod.o util_tvab_mod.o util_tvertical_geom_mod.o util_tveta_mod.o util_tvfe_mod.o util_txyb_mod.o util_txybder_mod.o util_yommp0_mod.o verder.o verdisint.o verint.o xrd_getoptions.o xrd_unix_env.o yemdim.o yemgeo.o yemgsl.o yemlap.o yemlbc_geo.o yemmp.o yoeaerc.o yoecmip.o yoeozoc.o yoeradghg.o yomcsgeom.o yomcst.o yomcver.o yomdim.o yomdimv.o yomdyn.o yomdyna.o yomgem.o yomgsgeom.o yomlap.o yomlddh.o yomleg.o yomlun.o yomlun_ifsaux.o yommp.o yommp0.o yomorog.o yomrip.o yomspddh.o yomsta.o yomtag.o yomvert.o
	$(FC) -o spcm.x spcm.o abor1.o geometry_mod.o intdyn_mod.o model_diagnostics_mod.o model_dynamics_mod.o model_general_conf_mod.o mxmaop.o mxptma.o mxture.o mxturs.o myrecvset.o mysendset.o reglatlon_field_mix.o set2pe.o sgemmx.o sigam_gp.o sigam_sp_openmp.o simplico.o sitnu_gp.o sitnu_sp_openmp.o spcimpfsolve.o spcm_simple.o spcsi.o trmtos.o trstom.o type_geometry.o type_model.o type_spgeom.o util_geometry_mod.o util_model_diagnostics_type_mod.o util_model_dynamics_type_mod.o util_model_general_conf_type_mod.o util_model_mod.o util_reglatlon_field_mod.o util_tcsgeom_blocked_mod.o util_tcsgeom_mod.o util_tcsgleg_mod.o util_tcst_mod.o util_tcty_mod.o util_tcver_mod.o util_tdim_mod.o util_tdimv_mod.o util_tdyn_mod.o util_tdyna_mod.o util_teaerc_macc_mod.o util_teaerc_tegen_mod.o util_tecmip_mod.o util_tedim_mod.o util_tegeo_mod.o util_tegsl_mod.o util_telbc_geo_mod.o util_temmp_mod.o util_teozoc_mod.o util_tgem_mod.o util_tgflt_mod.o util_tgmvt_mod.o util_tgsgeom_blocked_mod.o util_tgsgeom_mod.o util_thwind_mod.o util_tlap_mod.o util_tlddh_mod.o util_tlep_mod.o util_tmp_mod.o util_torog_blocked_mod.o util_torog_mod.o util_tpg_type_mod.o util_tradghg_mod.o util_trcp_mod.o util_trip_mod.o util_tspddh_mod.o util_tspgeom_mod.o util_tsta_mod.o util_ttnd_mod.o util_tvab_mod.o util_tvertical_geom_mod.o util_tveta_mod.o util_tvfe_mod.o util_txyb_mod.o util_txybder_mod.o util_yommp0_mod.o verder.o verdisint.o verint.o xrd_getoptions.o xrd_unix_env.o yemdim.o yemgeo.o yemgsl.o yemlap.o yemlbc_geo.o yemmp.o yoeaerc.o yoecmip.o yoeozoc.o yoeradghg.o yomcsgeom.o yomcst.o yomcver.o yomdim.o yomdimv.o yomdyn.o yomdyna.o yomgem.o yomgsgeom.o yomlap.o yomlddh.o yomleg.o yomlun.o yomlun_ifsaux.o yommp.o yommp0.o yomorog.o yomrip.o yomspddh.o yomsta.o yomtag.o yomvert.o $(LIBS)


subclean:
	\rm -f abor1.o mxmaop.o mxptma.o mxture.o mxturs.o myrecvset.o mysendset.o set2pe.o sgemmx.o sigam_gp.o sigam_sp_openmp.o simplico.o sitnu_gp.o sitnu_sp_openmp.o spcimpfsolve.o spcm.o spcm_simple.o spcsi.o trmtos.o trstom.o verder.o verdisint.o verint.o

fyppclean: 
	\rm -f sigam_gp.F90 sigam_sp_openmp.F90 sitnu_gp.F90 sitnu_sp_openmp.F90

clean: 
	\rm -f *.o *.xml *.a *.x *.mod *.optrpt 

tidy:
	\rm -f *.xml *.optrpt

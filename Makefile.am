bin_PROGRAMS = conf_check \
            conf_txt_to_bin \
            debug_rng \
            debug_su2 \
            debug_sun \
            debug_u1 \
            debug_twist \
            debug_agf_vs_gf \
            debug_agf_vs_delta \
            yang_mills_local_agf \
            yang_mills_local_pt_agf \
            yang_mills_conf_analysis \
            yang_mills_tuning_pt_mc \
            yang_mills_multilevel \
            yang_mills_multilevel_long \
            yang_mills_tracedef

timing = ./include/timing.h \
        ./lib/timing.c

align = ./include/memalign.h \
        ./lib/memalign.c

dSFMT = ./include/dSFMT-common.h \
        ./include/dSFMT.h \
        ./include/dSFMT-params11213.h \
        ./include/dSFMT-params1279.h \
        ./include/dSFMT-params132049.h \
        ./include/dSFMT-params19937.h \
        ./include/dSFMT-params216091.h \
        ./include/dSFMT-params2203.h \
        ./include/dSFMT-params4253.h \
        ./include/dSFMT-params44497.h \
        ./include/dSFMT-params521.h \
        ./include/dSFMT-params86243.h \
        ./include/dSFMT-params.h \
        ./lib/dSFMT.c

endian = ./include/endianness.h  \
         ./lib/endianness.c

gauge_group = ./include/gauge_group.h

gauge_conf = ./include/gauge_conf.h \
            ./lib/gauge_conf_def.c \
            ./lib/gauge_conf_meas.c \
            ./lib/gauge_conf_multilevel.c \
            ./lib/gauge_conf_upd.c \
            ./lib/gauge_conf_paral_temp.c \
            ./lib/gauge_conf_multicanonic.c

geometry = ./include/geometry.h \
           ./lib/geometry.c

gparam = ./include/gparam.h \
         ./lib/gparam.c

random = ./include/myrngparam.h \
         ./include/random.h \
         ./lib/random.c

su2 = ./include/su2.h \
      ./include/su2_upd.h

sun = ./include/sun.h \
      ./include/sun_upd.h

tens_prod = ./include/tens_prod.h

u1 = ./include/u1.h \
     ./include/u1_upd.h

conf_check_SOURCES=./src/conf_check.c \
                   $(timing) \
                   $(align) \
                   $(dSFMT) \
                   $(endian) \
                   $(gauge_group) \
                   $(gparam) \
                   $(random) \
                   $(su2) \
                   $(sun) \
                   $(u1)

conf_txt_to_bin_SOURCES=./src/conf_txt_to_bin.c \
                   $(timing) \
                   $(align) \
                   $(dSFMT) \
                   $(endian) \
                   $(gauge_group) \
                   $(gauge_conf) \
                   $(geometry) \
                   $(gparam) \
                   ./include/macro.h \
                   $(random) \
                   $(su2) \
                   $(sun) \
                   $(tens_prod) \
                   $(u1)

debug_rng_SOURCES=./src/debug_rng.c \
                  $(align) \
                  $(dSFMT) \
                  $(random)

debug_su2_SOURCES=./src/debug_su2.c \
                  $(timing) \
                  $(align) \
                  $(dSFMT) \
                  $(endian) \
                  $(gparam) \
                  ./include/macro.h \
                  $(random) \
                  $(su2)

debug_sun_SOURCES=./src/debug_sun.c \
                  $(timing) \
                  $(align) \
                  $(dSFMT) \
                  $(endian) \
                  $(gparam) \
                  ./include/macro.h \
                  $(random) \
                  $(su2) \
                  $(sun)

debug_u1_SOURCES=./src/debug_u1.c \
                  $(timing) \
                  $(align) \
                  $(dSFMT) \
                  $(endian) \
                  $(gparam) \
                  ./include/macro.h \
                  $(random) \
                  $(u1)

debug_twist_SOURCES=./src/debug_twist.c \
                          $(timing) \
                          $(align) \
                          $(dSFMT)  \
                          $(endian) \
                          $(gauge_group) \
                          $(gauge_conf) \
                          $(geometry) \
                          $(gparam) \
                          ./include/macro.h \
                          $(random) \
                          $(su2) \
                          $(sun) \
                          $(tens_prod) \
                          $(u1)

yang_mills_local_agf_SOURCES=./src/yang_mills_local_agf.c \
                          $(timing) \
                          $(align) \
                          $(dSFMT)  \
                          $(endian) \
                          $(gauge_group) \
                          $(gauge_conf) \
                          $(geometry) \
                          $(gparam) \
                          ./include/macro.h \
                          $(random) \
                          $(su2) \
                          $(sun) \
                          $(tens_prod) \
                          $(u1)

yang_mills_multilevel_SOURCES=./src/yang_mills_multilevel.c \
                             $(timing) \
                             $(align) \
                             $(dSFMT)  \
                             $(endian) \
                             $(gauge_group) \
                             $(gauge_conf) \
                             $(geometry) \
                             $(gparam) \
                             ./include/macro.h \
                             $(random) \
                             $(su2) \
                             $(sun) \
                             $(tens_prod) \
                             $(u1)

yang_mills_multilevel_long_SOURCES=./src/yang_mills_multilevel_long.c \
                                  $(timing) \
                                  $(align) \
                                  $(dSFMT)  \
                                  $(endian) \
                                  $(gauge_group) \
                                  $(gauge_conf) \
                                  $(geometry) \
                                  $(gparam) \
                                  ./include/macro.h \
                                  $(random) \
                                  $(su2) \
                                  $(sun) \
                                  $(tens_prod) \
                                  $(u1)

yang_mills_tracedef_SOURCES=./src/yang_mills_tracedef.c \
                            $(timing) \
                            $(align) \
                            $(dSFMT)  \
                            $(endian) \
                            $(gauge_group) \
                            $(gauge_conf) \
                            $(geometry) \
                            $(gparam) \
                            ./include/macro.h \
                            $(random) \
                            $(su2) \
                            $(sun) \
                            $(tens_prod) \
                            $(u1)

yang_mills_local_pt_agf_SOURCES=./src/yang_mills_local_pt_agf.c \
                          $(timing) \
                          $(align) \
                          $(dSFMT)  \
                          $(endian) \
                          $(gauge_group) \
                          $(gauge_conf) \
                          $(geometry) \
                          $(gparam) \
                          ./include/macro.h \
                          $(random) \
                          $(su2) \
                          $(sun) \
                          $(tens_prod) \
                          $(u1)

yang_mills_conf_analysis_SOURCES=./src/yang_mills_conf_analysis.c \
                       $(timing) \
                       $(align) \
                       $(dSFMT)  \
                       $(endian) \
                       $(gauge_group) \
                       $(gauge_conf) \
                       $(geometry) \
                       $(gparam) \
                       ./include/macro.h \
                       $(random) \
                       $(su2) \
                       $(sun) \
                       $(tens_prod) \
                       $(u1)

debug_agf_vs_gf_SOURCES=./src/debug_agf_vs_gf.c \
                          $(timing) \
                          $(align) \
                          $(dSFMT)  \
                          $(endian) \
                          $(gauge_group) \
                          $(gauge_conf) \
                          $(geometry) \
                          $(gparam) \
                          ./include/macro.h \
                          $(random) \
                          $(su2) \
                          $(sun) \
                          $(tens_prod) \
                          $(u1)

debug_agf_vs_delta_SOURCES=./src/debug_agf_vs_delta.c \
                          $(timing) \
                          $(align) \
                          $(dSFMT)  \
                          $(endian) \
                          $(gauge_group) \
                          $(gauge_conf) \
                          $(geometry) \
                          $(gparam) \
                          ./include/macro.h \
                          $(random) \
                          $(su2) \
                          $(sun) \
                          $(tens_prod) \
                          $(u1)

yang_mills_tuning_pt_mc_SOURCES=./src/yang_mills_tuning_pt_mc.c \
                          $(timing) \
                          $(align) \
                          $(dSFMT)  \
                          $(endian) \
                          $(gauge_group) \
                          $(gauge_conf) \
                          $(geometry) \
                          $(gparam) \
                          ./include/macro.h \
                          $(random) \
                          $(su2) \
                          $(sun) \
                          $(tens_prod) \
                          $(u1)

AM_CFLAGS = -std=c99 -Wall -Wextra -Werror -pedantic -Wconversion @OPENMP_CFLAGS@

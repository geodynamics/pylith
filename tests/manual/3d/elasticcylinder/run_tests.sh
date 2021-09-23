#!/bin/bash
pylith wedge_pres_elastic.cfg wedge_pres_elastic_b1_q1_hex.cfg > wedge_pres_elastic_b1_q1_hex.log 2>&1
pylith wedge_pres_elastic.cfg wedge_pres_elastic_b1_q1_tet.cfg > wedge_pres_elastic_b1_q1_tet.log 2>&1
pylith wedge_pres_elastic.cfg wedge_pres_elastic_b2_q2_hex.cfg > wedge_pres_elastic_b2_q2_hex.log 2>&1
pylith wedge_pres_elastic.cfg wedge_pres_elastic_b2_q2_tet.cfg > wedge_pres_elastic_b2_q2_tet.log 2>&1
pylith wedge_pres_elastic.cfg wedge_pres_elastic_b3_q3_hex.cfg > wedge_pres_elastic_b3_q3_hex.log 2>&1
pylith wedge_pres_elastic.cfg wedge_pres_elastic_b3_q3_tet.cfg > wedge_pres_elastic_b3_q3_tet.log 2>&1

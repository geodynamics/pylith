#!/bin/bash
pylith wedge_pres_elastic.cfg wedge_pres_elastic_b1_q1_quad.cfg > wedge_pres_elastic_b1_q1_quad.log 2>&1
pylith wedge_pres_elastic.cfg wedge_pres_elastic_b1_q1_tri.cfg > wedge_pres_elastic_b1_q1_tri.log 2>&1
pylith wedge_pres_elastic.cfg wedge_pres_elastic_b2_q2_quad.cfg > wedge_pres_elastic_b2_q2_quad.log 2>&1
pylith wedge_pres_elastic.cfg wedge_pres_elastic_b2_q2_tri.cfg > wedge_pres_elastic_b2_q2_tri.log 2>&1
pylith wedge_pres_elastic.cfg wedge_pres_elastic_b3_q3_quad.cfg > wedge_pres_elastic_b3_q3_quad.log 2>&1
pylith wedge_pres_elastic.cfg wedge_pres_elastic_b3_q3_tri.cfg > wedge_pres_elastic_b3_q3_tri.log 2>&1

#!/usr/bin/bash

export PYTHONPATH=$(dirname $(realpath ${BASH_SOURCE[0]}))
export OREKITDATA=$(locate orekit-data.zip)

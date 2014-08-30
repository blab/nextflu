#!/bin/bash
ntpdate pool.ntp.org
nohup Xvfb :99 -shmem -screen 0 1366x768x16 &
nohup x11vnc -display :99 -N -forever &
python augur/run.py

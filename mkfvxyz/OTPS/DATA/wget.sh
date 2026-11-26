#!/bin/bash

rsync -avz /mnt/share/Vol18T/nonoyama_tide/Kikai_OTPS/mkfvxyz/OTPS/DATA/ ./
#rsync -e ssh -avz minami@10.35.22.56:/mnt/share/Vol18T/nonoyama_tide/Kikai_OTPS/mkfvxyz/OTPS/DATA/ ./
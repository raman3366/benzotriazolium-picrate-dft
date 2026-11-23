#!/bin/bash
# R007 BTPC revision: one on-demand t2d-standard-32 VM for the new DFT runs (plan §5, decisions table D3 revised).
# Dry run by default; pass --yes to execute. Account, project and prices are recorded in the project notes (not deposited).
# Network: own VPC with SSH from this Pi's public address only (the project's default network allows SSH from anywhere).
set -euo pipefail
CFG=${CFG:-default}
ZONE=${ZONE:-us-west1-b}
REGION=${ZONE%-*}
NAME=${NAME:-r007-btpc-1}
MT=${MT:-t2d-standard-32}
DISK_GB=${DISK_GB:-200}
RUN_H=${RUN_H:-72}
KEY=${KEY:-$HOME/.ssh/r007_gcp_ed25519}
DO=echo; [ "${1:-}" = "--yes" ] && DO=
G="gcloud --configuration=$CFG --quiet"
PI_IP=$(curl -s --max-time 10 https://api.ipify.org)
[ -f "$KEY" ] || ssh-keygen -t ed25519 -N '' -C r007-btpc -f "$KEY" >/dev/null
echo "account=$CFG zone=$ZONE name=$NAME type=$MT disk=${DISK_GB}GB pd-balanced hard-stop=${RUN_H}h pi=$PI_IP/32"
$G compute networks describe r007-vpc >/dev/null 2>&1 || $DO $G compute networks create r007-vpc --subnet-mode=custom
$G compute networks subnets describe r007-sn-$REGION --region=$REGION >/dev/null 2>&1 || \
  $DO $G compute networks subnets create r007-sn-$REGION --network=r007-vpc --region=$REGION --range=10.97.0.0/24
$G compute firewall-rules describe r007-allow-ssh-pi >/dev/null 2>&1 || \
  $DO $G compute firewall-rules create r007-allow-ssh-pi --network=r007-vpc --direction=INGRESS --allow=tcp:22 \
     --source-ranges="$PI_IP/32" --target-tags=r007-vm
$DO $G compute instances create "$NAME" --zone="$ZONE" --machine-type="$MT" \
  --provisioning-model=STANDARD --max-run-duration="${RUN_H}h" --instance-termination-action=STOP \
  --network=r007-vpc --subnet=r007-sn-$REGION --tags=r007-vm --no-service-account --no-scopes \
  --image-family=ubuntu-2404-lts-amd64 --image-project=ubuntu-os-cloud \
  --boot-disk-size="${DISK_GB}GB" --boot-disk-type=pd-balanced \
  --labels=project=r007-btpc,owner=anantha \
  --metadata=ssh-keys="r007:$(cat "$KEY.pub")",block-project-ssh-keys=true

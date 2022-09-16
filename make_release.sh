#!/bin/sh
VER=`python3 version.py`
DCCVER="CPTAC-mzIdentML-v$VER"
./ulbin.sh
./ulsrc.sh
rm -f dist/cptacmzid-$VER.txt
( cd dist; md5sum cptacmzid-$VER.*.tgz > cptacmzid-$VER.md5 ; touch cptacmzid-$VER.txt )
if [ "$1" ]; then 
  for comment in "$@"; do 
    echo "* $comment" >> dist/cptacmzid-$VER.txt
  done
fi
gh release delete "$DCCVER" -y
git push --delete origin "refs/tags/$DCCVER"
git tag --delete "$DCCVER"
gh release create -F dist/cptacmzid-$VER.txt "$DCCVER" dist/cptacmzid-$VER.*.tgz dist/cptacmzid-$VER.md5
for a in dist/cptacmzid-$VER.*.tgz; do
  a1=`basename $a`
  rclone copyto $a cptac-s3:cptac-cdap.georgetown.edu/$a1
  b1=`echo $a1 | sed "s/-$VER//"`
  aws --profile cptac s3api put-object --bucket cptac-cdap.georgetown.edu --key "$b1" --website-redirect-location "/$a1"
done

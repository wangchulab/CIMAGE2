export CIMAGE_PATH=$(dirname $0)/..
for c in 2 3 4 5 6
do
  python $CIMAGE_PATH/python/match_paired_features.py $1 $c $2
done

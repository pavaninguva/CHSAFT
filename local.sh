# This is the postprocessing script

foldername=$(date +%Y%m%d_%H%M%S);
mkdir "$foldername";
cp params.yml ./"$foldername";
mv output* ./"$foldername";
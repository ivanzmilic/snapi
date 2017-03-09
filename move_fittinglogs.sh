
dest="$1"
suffix="$2"

for f in "$src"/*
do
       filename=$(basename "$f")

       filename="$dest/$filename-$suffix"

       echo "Copying $f to $filename"
       cp "$f" "$filename"
done

mv fitting_log.txt where/
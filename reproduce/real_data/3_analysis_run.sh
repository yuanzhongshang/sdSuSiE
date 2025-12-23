cd /net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8
awk 'NR > 3 && FNR <= 628 { print $2 }' /net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/trait.txt | while read value; do
if [ -e "$value/region.txt" ]; then
echo "value is: $value"
cp -r /net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/4079_irnt/3_analysis.sh /net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/$value/
cp -r /net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/4079_irnt/3_analysis.R /net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/$value/
cd /net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/$value/
sed -i "s/4079_irnt/$value/g" 3_analysis.sh
line_count=$(awk 'NR > 1 {count++} END {print count}' region.txt)
sed -i "s/1%/$line_count%/g" 3_analysis.sh
sbatch 3_analysis.sh
fi
cd /net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8
done

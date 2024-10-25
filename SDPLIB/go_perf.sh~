FLAMEGRAPH_DIR=/home/docker/FlameGraph

SDPxzs="\
arch0.dat-s.xz \
arch2.dat-s.xz \
arch4.dat-s.xz \
arch8.dat-s.xz \
control10.dat-s.xz \
control11.dat-s.xz \
control1.dat-s.xz \
control2.dat-s.xz \
control3.dat-s.xz \
control4.dat-s.xz \
control5.dat-s.xz \
control6.dat-s.xz \
control7.dat-s.xz \
control8.dat-s.xz \
control9.dat-s.xz \
equalG11.dat-s.xz \
equalG51.dat-s.xz \
gpp100.dat-s.xz \
gpp124-1.dat-s.xz \
gpp124-2.dat-s.xz \
gpp124-3.dat-s.xz \
gpp124-4.dat-s.xz \
gpp250-1.dat-s.xz \
gpp250-2.dat-s.xz \
gpp250-3.dat-s.xz \
gpp250-4.dat-s.xz \
gpp500-1.dat-s.xz \
gpp500-2.dat-s.xz \
gpp500-3.dat-s.xz \
gpp500-4.dat-s.xz \
hinf10.dat-s.xz \
hinf11.dat-s.xz \
hinf12.dat-s.xz \
hinf13.dat-s.xz \
hinf14.dat-s.xz \
hinf15.dat-s.xz \
hinf1.dat-s.xz \
hinf2.dat-s.xz \
hinf3.dat-s.xz \
hinf4.dat-s.xz \
hinf5.dat-s.xz \
hinf6.dat-s.xz \
hinf7.dat-s.xz \
hinf8.dat-s.xz \
hinf9.dat-s.xz \
infd1.dat-s.xz \
infd2.dat-s.xz \
infp1.dat-s.xz \
infp2.dat-s.xz \
maxG11.dat-s.xz \
maxG32.dat-s.xz \
maxG51.dat-s.xz \
maxG55.dat-s.xz \
maxG60.dat-s.xz \
mcp100.dat-s.xz \
mcp124-1.dat-s.xz \
mcp124-2.dat-s.xz \
mcp124-3.dat-s.xz \
mcp124-4.dat-s.xz \
mcp250-1.dat-s.xz \
mcp250-2.dat-s.xz \
mcp250-3.dat-s.xz \
mcp250-4.dat-s.xz \
mcp500-1.dat-s.xz \
mcp500-2.dat-s.xz \
mcp500-3.dat-s.xz \
mcp500-4.dat-s.xz \
qap10.dat-s.xz \
qap5.dat-s.xz \
qap6.dat-s.xz \
qap7.dat-s.xz \
qap8.dat-s.xz \
qap9.dat-s.xz \
qpG11.dat-s.xz \
qpG51.dat-s.xz \
ss30.dat-s.xz \
theta1.dat-s.xz \
theta2.dat-s.xz \
theta3.dat-s.xz \
theta4.dat-s.xz \
theta5.dat-s.xz \
theta6.dat-s.xz \
thetaG11.dat-s.xz \
thetaG51.dat-s.xz \
truss1.dat-s.xz \
truss2.dat-s.xz \
truss3.dat-s.xz \
truss4.dat-s.xz \
truss5.dat-s.xz \
truss6.dat-s.xz \
truss7.dat-s.xz \
truss8.dat-s.xz"

for SDPxz in $SDPxzs; do
    SDP="${SDPxz%.xz}"
    rm -f $SDP $SDPxz
    cp data/$SDPxz .
    unxz $SDPxz
    start_time=$(date +%s%N)
    sudo perf record -o perf.data_${SDP%.*} -g -- ../sdpa_gmp -ds $SDP -o ${SDP%.*}.result
    end_time=$(date +%s%N)
    duration=$(( (end_time - start_time) / 1000000 ))
    echo "Execution time: $duration ms"

    sudo perf script -i perf.data_${SDP%.*} > out_${SDP%.*}.perf
    cat out_${SDP%.*}.perf | $FLAMEGRAPH_DIR/stackcollapse-perf.pl | $FLAMEGRAPH_DIR/flamegraph.pl > flamegraph_${SDP%.*}.svg
    echo "Flamegraph for $exe generated."

    sudo perf annotate -i perf.data_${SDP%.*} > annotation_${SDP%.*}.txt
    echo "Annotation for $exe saved."
    sudo perf report -i perf.data_${SDP%.*} --stdio > report_${SDP%.*}.txt
    echo "Perf report for $exe generated and saved to report_${SDP%.*}.txt."
    echo
done

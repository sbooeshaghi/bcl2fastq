#!/bin/bash

#set -e
#set -x

# number of cycles
cycles_count=200
# first index cycle
cycle_index_start=92
# last index cycle
cycle_index_end=108
# number of clusters
clusters_count=10000
# tiles aggregation (1: nova, 0: hiseq/miseq)
aggregate_tiles=1

# generates predictable data with following properties:
# - clusters are assigned to samples A,C,G,T periodically and in this order 
# - all clusters for samples C and G have filter flag set to true 
# - all clusters for sample T have filter flag set to false 
# - filter flag for clusters of sample A is alternating between true and false periodically
# - there are 10000 clusters on 9 tiles 
# - tiles contain following number of clusters: 3,9,27,81,etc. (since 10000 is not divisible by 3, 9th tile contains only 160 clusters)
# - sample sheet, runinfo XML and config XML (last one only for hiseq/miseq) need to be created by hand (examples below are compatible with above described data)

# example sample sheet (<runfolder>/SampleSheet.csv)
#[Header],,,,
#Investigator Name,Isabelle,,,
#Project Name,Nova,,,
#Experiment Name,Orbital death ray research volume LXXIV,,,
#Date,5/27/2025,,,
#Workflow,GenerateFASTQ,,,
#,,,,
#[Settings],,,,
#MaskAdapter,CGCGTATACGCGTATA,,,
#TrimAdapter,GCGCATATGCGCATAT,,,
#,,,,
#[Data],,,,
#SampleID,SampleName,index,index2
#AA,AA,AAAAAAAA,AAAAAAAA
#CC,CC,CCCCCCCC,CCCCCCCC
#GG,GG,GGGGGGGG,GGGGGGGG
#TT,TT,TTTTTTTT,TTTTTTTT
#XX,XX,ACGTACGT,ACGTACGT

# example runinfo XML (<runfolder>/RunInfo.xml)
#<?xml version="1.0"?>
#<RunInfo xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" Version="2">
#  <Run Id="Test" Number="42">
#    <Flowcell>Test</Flowcell>
#    <Instrument>CSSIM</Instrument>
#    <Date>8/20/2013 3:29:17 PM</Date>
#    <Reads>
#      <Read FirstCycle="1" LastCycle="91" IsIndexedRead="N" />
#      <Read FirstCycle="92" LastCycle="99" IsIndexedRead="Y" />
#      <Read FirstCycle="100" LastCycle="107" IsIndexedRead="Y" />
#      <Read FirstCycle="108" LastCycle="200" IsIndexedRead="N" />
#    </Reads>
#    <FlowcellLayout LaneCount="1" SurfaceCount="2" SwathCount="3" TileCount="6" SectionPerLane="3" LanePerSection="2" />
#  </Run>
#</RunInfo>

# example config XML (<runfolder>/Data/Intensities/BaseCalls/config.xml); hiseq/miseq only
#<?xml version="1.0" encoding="utf-8"?>
#<BaseCallAnalysis xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
#  <Run Name="BaseCalls">
#    <TileSelection>
#      <Lane Index="1">
#        <Tile>1</Tile>
#        <Tile>2</Tile>
#        <Tile>3</Tile>
#        <Tile>4</Tile>
#        <Tile>5</Tile>
#        <Tile>6</Tile>
#        <Tile>7</Tile>
#        <Tile>8</Tile>
#        <Tile>9</Tile>
#      </Lane>
#    </TileSelection>
#  </Run>
#</BaseCallAnalysis>


mkdir -p ./Data/Intensities/BaseCalls/L001/
mkdir -p ./Data/Intensities/L001/

tmp_bcl_filename='./tmp_bcl'
bci_filename='./Data/Intensities/BaseCalls/L001/s_1.bci'
rm -f "$bci_filename"
locs_filename='./Data/Intensities/L001/s_1.locs'
filter_filename='./Data/Intensities/BaseCalls/L001/s_1.filter'

for cycle_no in $(seq 1 $cycles_count)
do
  echo "Generating cycle: $cycle_no"

  if [ $aggregate_tiles -eq 1 ]
  then
    bcl_filename=`printf './Data/Intensities/BaseCalls/L001/%04u.bcl.bgzf' $cycle_no`
    printf '0: %.8x' $clusters_count | sed -E 's/0: (..)(..)(..)(..)/0: \4\3\2\1/' | xxd -r -g0 > "$tmp_bcl_filename"
    gzip "$tmp_bcl_filename"
    cat "$tmp_bcl_filename.gz" > "$bcl_filename"
    rm "$tmp_bcl_filename.gz"

    if [ $cycle_no -eq 1 ]
    then
      printf '0: 010000000000803f' | xxd -r -g0 > "$locs_filename"
      printf '0: %.8x' $clusters_count | sed -E 's/0: (..)(..)(..)(..)/0: \4\3\2\1/' | xxd -r -g0 >> "$locs_filename"

      printf '0: 0000000003000000' | xxd -r -g0 > "$filter_filename"
      printf '0: %.8x' $clusters_count | sed -E 's/0: (..)(..)(..)(..)/0: \4\3\2\1/' | xxd -r -g0 >> "$filter_filename"
    fi
  else
    bcl_dir=`printf './Data/Intensities/BaseCalls/L001/C%d.1/' $cycle_no`
    mkdir -p "$bcl_dir"
  fi


  tile_idx=0
  cluster_idx=0
  clusters_on_current_tile=3
  while [ $cluster_idx -lt $clusters_count ]
  do
    tile_cluster_idx=0
    while [[ $cluster_idx -lt $clusters_count ]] && [[ $tile_cluster_idx -lt $clusters_on_current_tile ]]
    do
      if [[ $cycle_no -ge $cycle_index_start ]] && [[ $cycle_no -lt $cycle_index_end ]]
      then
        base=$(($cluster_idx % 4))
        quality=0x3F
      else
        base=$(($(($cluster_idx + $cycle_no - 1)) % 0x04))
        quality=$(($(($cluster_idx + $cycle_no - 1)) % 0x3F))
      fi
      printf '0: %.2x' $(($(($quality << 2)) | $base)) | xxd -r -g0 >> "$tmp_bcl_filename"

      if [ $cycle_no -eq 1 ]
      then
        printf '0: cdcc8c3f9a99993f' | xxd -r -g0 >> "$locs_filename"
        filter=1;
        if [[ $(($cluster_idx % 8)) -eq 0 ]] || [[ $(($cluster_idx % 4)) -eq 3 ]]
        then
            filter=0;
        fi
        printf '0: %.2x' $filter | xxd -r -g0 >> "$filter_filename"
        #printf '0: 01' | xxd -r -g0 >> "$filter_filename"
      fi

      tile_cluster_idx=$(($tile_cluster_idx + 1))
      cluster_idx=$(($cluster_idx + 1))
    done
    if [ $aggregate_tiles -eq 1 ]
    then
        gzip "$tmp_bcl_filename"
        cat "$tmp_bcl_filename.gz" >> "$bcl_filename"
        rm "$tmp_bcl_filename.gz"
    else
      tile_number=$(($tile_idx + 1))
      bcl_filename=`printf '%s/s_1_%d.bcl.gz' $bcl_dir $tile_number`
      printf '0: %.8x' $tile_cluster_idx | sed -E 's/0: (..)(..)(..)(..)/0: \4\3\2\1/' | xxd -r -g0 > "$tmp_bcl_filename.h"
      cat "$tmp_bcl_filename" >> "$tmp_bcl_filename.h"
      rm "$tmp_bcl_filename"
      gzip "$tmp_bcl_filename.h"
      cat "$tmp_bcl_filename.h.gz" > "$bcl_filename"
      rm "$tmp_bcl_filename.h.gz"

      if [ $cycle_no -eq 1 ]
      then
        correct_filter_filename=`printf "./Data/Intensities/BaseCalls/L001/s_1_%d.filter" $tile_number`
        printf '0: %.8x' $tile_cluster_idx | sed -E 's/0: (..)(..)(..)(..)/0: \4\3\2\1/' | xxd -r -g0 > "$correct_filter_filename"
        cat "$filter_filename" >> "$correct_filter_filename"
        rm "$filter_filename"

        correct_locs_filename=`printf "./Data/Intensities/L001/s_1_%d.locs" $tile_number`
        printf '0: %.8x' $tile_cluster_idx | sed -E 's/0: (..)(..)(..)(..)/0: \4\3\2\1/' | xxd -r -g0 > "$correct_locs_filename"
        cat "$locs_filename" >> "$correct_locs_filename"
        rm "$locs_filename"
      fi
    fi

    tile_idx=$(($tile_idx + 1))
    if [ $cycle_no -eq 1 ]
    then
      printf '0: %.8x' $tile_idx | sed -E 's/0: (..)(..)(..)(..)/0: \4\3\2\1/' | xxd -r -g0 >> "$bci_filename"
      printf '0: %.8x' $tile_cluster_idx | sed -E 's/0: (..)(..)(..)(..)/0: \4\3\2\1/' | xxd -r -g0 >> "$bci_filename"
    fi

    clusters_on_current_tile=$(($clusters_on_current_tile * 3))
  done
done



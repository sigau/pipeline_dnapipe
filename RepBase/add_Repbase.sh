#! /bin/bash
cd /mnt2/new_lib/
cp /opt/RepeatMasker/Libraries/* /mnt2/new_lib/
cp /mnt2/tmp_repbase/* /mnt2/new_lib/ 
addRepBase.pl -libdir /mnt2/new_lib

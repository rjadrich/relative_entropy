../multi_g_rdf -b 200 -e 798 -- -bin 0.005 << 'EOF' 
2
2
'EOF'

rdf_exit=$?

exit $rdf_exit

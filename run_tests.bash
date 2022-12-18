export EPICS_CA_MAX_ARRAY_BYTES=1000000000000
python tests/start_server.py &
sleep 5
pytest tests

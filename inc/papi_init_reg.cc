#ifdef PAPI
	int event_set;
	long_long values[4];
	int retval; 
	event_set = PAPI_NULL;
	values [0] = values[1] = values [2] = values[3] = 0;

	retval = PAPI_library_init(PAPI_VER_CURRENT);

	if(retval != PAPI_VER_CURRENT && retval > 0){
		cout<<retval << " "<< PAPI_VER_CURRENT<<endl;
	}

	string event_names[4];

	event_names[0] = event_names[1] = event_names[2] = event_names[3] = "event";


	if ((retval = PAPI_create_eventset(&event_set)) != PAPI_OK) {
		fprintf(stderr, "PAPI error %d: %s\n",retval,           PAPI_strerror(retval));
		exit(1);
	}

	//PAPI_ERROR_CHECK(PAPI_add_event(event_set, PAPI_L1_DCM));
	PAPI_ERROR_CHECK(PAPI_add_event(event_set, PAPI_LD_INS));
	event_names[0] += "PAPI_LD_INS";
	PAPI_ERROR_CHECK(PAPI_add_event(event_set, PAPI_SR_INS));
	event_names[1] += "PAPI_SR_INS";
	PAPI_ERROR_CHECK(PAPI_add_event(event_set, PAPI_L1_DCM));
	event_names[2] += "PAPI_L1_DCM";

	PAPI_ERROR_CHECK(PAPI_start(event_set));
#endif
#ifdef PAPI
	PAPI_ERROR_CHECK(PAPI_stop(event_set, values));
	std::cout << event_names[0]<<",  " << values[0] << ", " <<  event_names[1]<<",  "  << values[1]<<", " << event_names[2]<<",  "  << values[2]<<", " << event_names[3]<<",  "  << values[3] << "\n";
#endif
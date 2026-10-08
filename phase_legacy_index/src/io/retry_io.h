/*******************************************************************************
 * Copyright (C) 2022-2023 Simone Rubinacci
 * Copyright (C) 2022-2023 Olivier Delaneau
 *
 * MIT Licence
 ******************************************************************************/

#ifndef _RETRY_IO_H
#define _RETRY_IO_H

#include <string>
#include <chrono>
#include <thread>

#include <utils/otools.h>

struct attempt_result
{
	bool ok = false;
	bool fatal = false;
	std::string err;
};

template <typename Fn>
void retry_with_backoff(const std::string& what, int n_retry, std::chrono::seconds base_delay, Fn&& attempt)
{
	std::chrono::seconds delay = base_delay;
	for (int i = 0; i < n_retry; ++i)
	{
		const attempt_result r = attempt();
		if (r.ok) return;
		if (r.fatal) vrb.error("Non-retryable error while " + what + ": " + r.err);

		const int attempts_left = n_retry - i - 1;
		if (attempts_left > 0)
		{
			vrb.warning("Error while " + what + ": " + r.err + ". " + std::to_string(attempts_left) + " attempt(s) remaining");
			std::this_thread::sleep_for(delay);
			delay *= 2;
		}
		else vrb.error("Max number of retries attempted while " + what + ". Last error: " + r.err);
	}
}

#endif

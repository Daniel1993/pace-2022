#ifndef AUX_H_GUARD
#define AUX_H_GUARD

#ifndef NDEBUG
#define DEB_PRINTF(...) printf(__VA_ARGS__)
#else /* DEBUG */
#define DEB_PRINTF(...) /* empty */
#endif /* NDEBUG */

#define TIMER_T struct timespec
#define TIMER_READ(_timer) clock_gettime(CLOCK_MONOTONIC_RAW, (struct timespec*)&(_timer));

#define CLOCK_DIFF_MS(startClock, endClock) \
  ((double)(endClock.tv_sec-startClock.tv_sec) * 1000.0f + \
  (double)(endClock.tv_nsec-startClock.tv_nsec) / 1.0e6f)

#include <errno.h>
#include <stdio.h>
#include <signal.h>
#include <string.h>
#include <unistd.h>

extern volatile int sigterm_called; // set to one after sigterm

static void sigterm_handler(int sig, siginfo_t *si, void *unused)
{
  sigterm_called = 1;
}

static void handle_sigterm()
{
  struct sigaction sa;
  memset(&sa, 0, sizeof(sa));
  sigemptyset(&sa.sa_mask);
  sa.sa_sigaction = sigterm_handler;
  if (-1 == sigaction(/* soft kill */SIGALRM, &sa, NULL))
    perror("sigaction(SIGALRM)");
  if (-1 == sigaction(/* soft kill */SIGTERM, &sa, NULL))
    perror("sigaction(SIGTERM)");
  if (-1 == sigaction(/* ctrl+c */SIGINT, &sa, NULL))
    perror("sigaction(SIGINT)");
}

#endif /* AUX_H_GUARD */
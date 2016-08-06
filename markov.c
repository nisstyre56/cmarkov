#include <time.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include "error.h"
#include "maa.h"
#include "tokenize.h"
#include "markov.h"

#define LEN 50

static char quote = '\'';

static inline void
initialize_neighbours(const char *str,
                     hsh_HashTable graph_table) {
  /* Initialize the table of neighbours corresponding to some string in the markov graph */
  assert(!hsh_retrieve(graph_table, str));
  neighbours_t *new_neighbours = xmalloc(sizeof (neighbours_t));
  CHECK(new_neighbours);
  new_neighbours->neighbours = hsh_create(NULL, NULL);
  new_neighbours->number = 0;
  new_neighbours->unique_num = 0;
  hsh_insert(graph_table, str, new_neighbours);
  return;
}

static inline unsigned long
numberof_keys(graph_t graph) {
  /* Get the number of unique keys in the graph */
  hsh_Stats stats = hsh_get_stats(graph.cache);
  unsigned long num = stats->entries;
  xfree(stats);
  return num;
}

static inline unsigned long
numberof_transitionable(graph_t graph) {
  /* Get the number of keys with >0 neighbours */
  /* Only call after graph has been converted */
  unsigned long num = 0;
  void *p, *key;
  markov_trans_t *val;
  HSH_ITERATE(graph.graph, p, key, val) {
    if (val->number > 0) {
      num++;
    }
  }
  return num;
}

static const char*
get_ngram(const char* str,
          graph_t graph) {
  /* Try to get a string from the cache.
   * If it's not already cached, allocate the memory for it
   * then return the freshly cached string
   */
  hsh_HashTable cache = graph.cache;
  hsh_HashTable graph_table = graph.graph;
  const char *exists = hsh_retrieve(cache, str);
  if (exists) {
    return exists;
  }
  else {
    /* Add it to the cache and return it */
    size_t gram_size = strlen(str) + 1;
    char *new_str = xmalloc(gram_size);
    CHECK(str);
    snprintf(new_str, gram_size, "%s", str);
    hsh_insert(cache, new_str, new_str);
    initialize_neighbours(new_str, graph_table);
    return new_str;
  }
}

static inline void
insert_neighbour(const char *left,
                 const char *neighbour,
                 graph_t graph) {
  /* Insert a neighbour into the table of neighbours for a given key */
  neighbours_t *neighbours = (neighbours_t *)hsh_retrieve(graph.graph, left);
  neighbours->number++;
  hsh_HashTable neighbours_table = neighbours->neighbours;
  CHECK(neighbours->neighbours);
  if (hsh_retrieve(neighbours_table, neighbour)) {
    return;
  }
  neighbours->unique_num++;
  const char *new_neighbour = get_ngram(neighbour, graph);
  CHECK(new_neighbour);
  uint32_t *count = xmalloc(sizeof (uint32_t));
  CHECK(count);
  *count = 0;
  hsh_insert(neighbours_table, new_neighbour, count);
}

static inline void
increment_neighbour(const char *left,
                    const char *neighbour,
                    graph_t graph) {
  /* Increment the frequency of a given bi-gram.
   * bi-gram does not necessarily mean a specific thing
   * it could be pairs of words, pairs of letters, sequences of n letters, and so on
   */
  neighbours_t *neighbours = (neighbours_t *)hsh_retrieve(graph.graph, left);
  hsh_HashTable neighbours_hash = neighbours->neighbours;
  CHECK(neighbours_hash);
  uint32_t *count = (uint32_t *)hsh_retrieve(neighbours_hash, neighbour);
  CHECK(count);
  (*count)++;
}

static inline neighbours_t*
get_neighbours(graph_t graph,
               char *gram) {
  /* Simply return the table of neighbours corresponding to a given string */
  neighbours_t *neighbours;
  neighbours = (neighbours_t *)hsh_retrieve(graph.graph, gram);
  assert(neighbours);
  return neighbours;
}

static inline markov_trans_t*
get_prob_neighbours(graph_t graph,
                    char *gram) {
  /* Return the converted probability transitions */
  markov_trans_t *neighbours;
  neighbours = (markov_trans_t *)hsh_retrieve(graph.graph, gram);
  assert(neighbours);
  return neighbours;
}

static inline void
convert_neighbours(graph_t graph,
                   char *gram) {
  neighbours_t *neighbours = get_neighbours(graph, gram);

  markov_trans_t *result = xmalloc(sizeof (markov_trans_t));
  CHECK(result);
  size_t nb_size = neighbours->number;
  hsh_HashTable neighbours_hash = neighbours->neighbours;

  void *key;
  uint32_t *frequency;
  void *p;
  uint32_t index = 0;
  probability_t transition;
  probability_t *neighbour_array = xcalloc(sizeof (probability_t), nb_size);
  CHECK(neighbour_array);
  HSH_ITERATE(neighbours_hash, p, key, frequency) {
    transition.frequent.frequency = *frequency;
    xfree(frequency);
    transition.frequent.token = key;
    neighbour_array[index] = transition;
    index++;
  }
 float lower = 0.0;
  probability_t current;
  for (uint32_t i = 0; i < neighbours->unique_num; i++) {
    current.frequent = neighbour_array[i].frequent;
    neighbour_array[i].bucket.token = current.frequent.token;
    neighbour_array[i].bucket.lower = lower;
    neighbour_array[i].bucket.upper = lower + ((float)neighbour_array[i].frequent.frequency) /
                                              (neighbours->number);
    lower = neighbour_array[i].bucket.upper;
  }
  result->transitions = neighbour_array;
  result->number = neighbours->unique_num;
  hsh_delete(graph.graph, gram);
  hsh_insert(graph.graph, gram, result);
  hsh_destroy(neighbours->neighbours);
  xfree(neighbours);
}

static inline void
convert_all_neighbours(graph_t graph) {
  void *p, *key;
  char *current_key;
  unsigned long num_keys = numberof_keys(graph);
  if (num_keys == 0) {
    return;
  }
  stk_Stack keys = stk_create();

  /* iterate over all keys K, in hash table T */
  HSH_ITERATE_KEYS(graph.graph, p, key) {
    stk_push(keys, key);
  }

  for (uint32_t i = 0; i < num_keys; i++) {
    current_key = (char *)stk_pop(keys);
    convert_neighbours(graph, current_key);
  }
  stk_destroy(keys);
}

static inline void
relate_bigram(const char *a,
              const char *b,
              graph_t graph) {
  /* Update the graph with the information that b follows a */
  const char* str = get_ngram(a, graph);
  insert_neighbour(str, b, graph);
  increment_neighbour(str, b, graph);
}

static int
transition_cmp(const void *keyval,
               const void *datum) {
  float chosen_number = *((float *)keyval);
  probability_t *transition = (probability_t *)datum;
  float lower = transition->bucket.lower;
  float upper = transition->bucket.upper;
  if ((chosen_number >= lower) &&
      (chosen_number <= upper)) {
    return 0;
  }
  else if (chosen_number < lower) {
    return -1;
  }
  else {
    return 1;
  }
}

static inline char*
pick_random_transition(unique_keys_t unique_neighbours) {
  size_t num = unique_neighbours.number;
  char **keys = unique_neighbours.keys;
  size_t selection = (size_t)floor(drand48() * (num - 1));
  return keys[selection];
}


static inline char*
next_ngram(graph_t graph,
           char *start,
           unique_keys_t unique_neighbours) {
  markov_trans_t *transitions = get_prob_neighbours(graph, start);
  if (transitions->number == 0) {
    return pick_random_transition(unique_neighbours);
  }
  probability_t *buckets = transitions->transitions;
  size_t bucket_size = transitions->number;
  float chosen = (float)drand48();
  probability_t *result = bsearch(&chosen,
                                  buckets,
                                  bucket_size,
                                  sizeof (probability_t),
                                  transition_cmp);
  return ((char *)result->bucket.token);
}

lst_List
generate_strings(markov_chain_t markov_chain,
                 char *start,
                 uint32_t n) {
  unique_keys_t unique_neighbours = markov_chain.unique;
  graph_t graph = markov_chain.graph;
  lst_List result = lst_create();
  char *current = start;
  for (uint32_t i = 0; i < n; i++) {
    lst_append(result, current);
    current = next_ngram(graph, current, unique_neighbours);
  }
  return result;
}

static inline unique_keys_t
get_all_keys(graph_t graph) {
  /* Gets all unique keys with neighbours */
  /* Should only be called after graph generation */
  unsigned long number = numberof_transitionable(graph);
  char **keys = xcalloc(sizeof (char *), number);
  CHECK(keys);
  void *p, *key;
  unique_keys_t result;
  markov_trans_t *val;
  uint32_t i = 0;
  HSH_ITERATE(graph.graph, p, key, val) {
    if (val->number > 0) {
      keys[i] = key;
      i++;
    }
  }
  result.keys = keys;
  result.number = i;
  return result;
}

static inline graph_t
make_graph(void) {
  /* Make an initial empty graph */
  graph_t result;
  result.cache = hsh_create(NULL, NULL);
  result.graph = hsh_create(NULL, NULL);
  return result;
}

static inline void
release_converted_graph(graph_t graph) {
  void *p, *key;
  markov_trans_t *datum;
  /* iterate over all keys K, in hash table and xfree them*/
  HSH_ITERATE(graph.graph, p, key, datum) {
    xfree(datum->transitions);
    xfree(datum);
    xfree(key);
  }
  hsh_destroy(graph.cache);
  hsh_destroy(graph.graph);
}

markov_chain_t
build_markov_chain(token_stream tokens) {
  markov_chain_t result;
  graph_t graph = make_graph();
  token_t current;
  token_t next;
  while (tokens.length > 1) {
    current = peek_token(&tokens);
    pop_token(&tokens);
    next = peek_token(&tokens);
    relate_bigram(token_to_string(next), token_to_string(current), graph);
  }
  convert_all_neighbours(graph);
  result.graph = graph;
  result.unique = get_all_keys(graph);
  return result;
}

char *
token_to_string(token_t token) {
    switch (token.token_type) {
      case WORD:
        return (char*)token.token.word;
        break;
      case INTEGER:
        return (char*)token.token.integer;
        break;
      case FLOATING:
        return (char*)token.token.floating;
        break;
      case QUOTE:
        return &quote;
        break;
      case PAREN:
        return (char*)token.token.parenthesis;
        break;
      case EMPTY:
        printf("should not be here\n");
        exit(EXIT_FAILURE);
        break;
      case STRING:
        return (char*)token.token.string;
        break;
      default:
        printf("oops, there was an unknown token, check valgrind or gdb\n");
        exit(EXIT_FAILURE);
    }
}

void
release_markov_chain(markov_chain_t chain) {
  release_converted_graph(chain.graph);
  xfree(chain.unique.keys);
  return;
}

int
main (void) {
  void *test_input = xmalloc(555000);
  size_t nbytes = read(STDIN_FILENO, test_input, 555000);

  if (nbytes == 0) {
    exit(EXIT_FAILURE);
  }
  token_stream test_bigrams_stack = tokenize(test_input, 0, nbytes);
  markov_chain_t chain = build_markov_chain(test_bigrams_stack);
  srand48(time(NULL));
  lst_List test = generate_strings(chain, token_to_string(peek_token(&test_bigrams_stack)), LEN);
  lst_pop(test);
  for (uint32_t i = 0; i < LEN-1; i++) {
    printf("%s ", (char *)lst_pop(test));
  }
  printf("\n");
  lst_destroy(test);
  _lst_shutdown();
  release_markov_chain(chain);
  xfree(test_input);
  release_tokens(&test_bigrams_stack);
  return EXIT_SUCCESS;
}

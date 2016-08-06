typedef
  struct {
    hsh_HashTable cache;
    hsh_HashTable graph;
  }
  graph_t;

typedef
  struct {
    size_t number;
    char **keys;
  }
  unique_keys_t;

typedef
  struct {
    graph_t graph;
    unique_keys_t unique;
  }
  markov_chain_t;


typedef
  struct {
    hsh_HashTable neighbours;
    size_t number;
    size_t unique_num;
  }
  neighbours_t;

/*
 * Transition types for various reasons
 */

typedef
  struct {
    float upper;
    float lower;
    const char *token;
  }
  bucket_t;

typedef
  struct {
    uint32_t frequency;
    const char *token;
  }
  transition_t;


typedef
  union {
    transition_t frequent;
    bucket_t bucket;
  }
  probability_t;

typedef
  struct {
    size_t number;
    probability_t *transitions;
  }
  markov_trans_t;


markov_chain_t
build_markov_chain(token_stream);

char *
token_to_string(token_t);

void
release_markov_chain(markov_chain_t);

lst_List
generate_strings(markov_chain_t,
                 char *,
                 uint32_t);

#define STACK_SIZE 4096
#define GROWTH_FACTOR 2
#define MAX_STRING_SIZE 30

typedef char* source_t;

typedef enum {
  WORD = 1,
  INTEGER = 2,
  FLOATING = 3,
  QUOTE = 4,
  WSPACE = 5,
  PAREN = 6 ,
  EMPTY = 7,
  STRING = 8
} tok_t;

typedef union {
    const char *word;
    const char *integer;
    const char *floating;
    const char *parenthesis;
    const char *string;
    bool quote;
    bool null_token;
} token_val_t;

typedef struct {
  tok_t token_type;
  token_val_t token;
} token_t;

typedef struct {
  size_t length; /* Number of current elements */
  size_t max_length; /* Maximum length of the stack */
  token_t *tokens;
  hsh_HashTable memo;
} token_stream;

bool
push_token(token_stream*, token_t);

bool
pop_token(token_stream*);

token_t
peek_token(token_stream*);

token_stream
tokenize(source_t, uint32_t, const uint32_t);

bool
release_tokens(token_stream*);

#ifndef TOK_LIB
static uint32_t
match_int(source_t, uint32_t, const uint32_t);

static uint32_t
match_float(source_t, uint32_t, const uint32_t);

static uint32_t
match_word(source_t, uint32_t, const uint32_t);

static uint32_t
match_string(source_t, uint32_t, const uint32_t);
#endif

int
free_token(const void *,
           const void *);
token_t
testfunc(void);

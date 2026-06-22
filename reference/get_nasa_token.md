# Get a NASA Earthdata bearer token

Fetch a NASA Earthdata bearer token from using the Earthdata API. If
none exist, this will create one, or if one already exists it will fetch
that one instead.

## Usage

``` r
get_nasa_token(username, password)
```

## Arguments

- username:

  character. NASA Earthdata username

- password:

  character. NASA Earthdata password

## Value

character

## Author

Simon E. H. Smart <simon.smart@cantab.net>

from config import CLIENT_DIST_PATH
#from BeautifulSoup import BeautifulSoup
from bs4 import BeautifulSoup
from shutil import copyfile
from argparse import ArgumentParser, REMAINDER
from re import sub, IGNORECASE

from codecs import encode

def main():

    args = parse_arguments()

    for arg in args.args:
        if args.backup:
            # save a backup of original just in case
            backup_name = arg + '.backup'
            copyfile(arg, backup_name)
            print ('Saved backup'), backup_name
        if args.remove:
            remove_script(arg)
        elif args.removeall:
            remove_all_scripts(arg)
        elif args.quick:
            add_preview_script_quick(arg)
        else:
            add_preview_script(arg)

        print ('Processed'), arg

    print ('Finished')

def get_script():
    # note that vendor folder only exists after successful build with gulp
    preview_script_src = '/'.join([CLIENT_DIST_PATH, 'vendor', 'require.js'])
    print(preview_script_src)
    soup = BeautifulSoup()
    script = soup.new_tag('script')
    script['src'] = preview_script_src
    script['data-main'] = '/'.join([CLIENT_DIST_PATH, 'main'])
    script['type'] = 'text/javascript'

    return script

def remove_script(file_name):
    script = str(get_script())
    raw = get_file_contents(file_name)
    processed = sub(script, '', raw, flags=IGNORECASE)
    write_new(file_name, processed)

def remove_all_scripts(file_name):
    raw = get_file_contents(file_name)
    processed = sub('<script.*</script>', '', raw, flags=IGNORECASE)
    write_new(file_name, processed)

def add_preview_script_quick(file_name):
    script = str(get_script())
    raw = get_file_contents(file_name)

    if raw.find(script) == -1:
        processed = sub('</head>', script + '</head>', raw, flags=IGNORECASE)
        write_new(file_name, processed)
    else:
        print ('Script already found in'), file_name

def add_preview_script(file_name):
    markup = get_file_contents(file_name)
    soup = BeautifulSoup(markup, 'html5lib')

    head = soup.head
    new_script = get_script()
    scripts = head.find_all('script', src = new_script.src)

    # check if script has already been added before adding script
    if not scripts:
        head.append(new_script)
        write_new(file_name, str(soup))
    else:
        print ('Script already found in'), file_name

def get_file_contents(file_name):
    with open(file_name, 'rb') as raw_file:
        return raw_file.read()

def write_new(file_name, contents):
    with open(file_name, 'w', encoding='utf-8') as f:
        f.write(contents)

def parse_arguments():
    parser = ArgumentParser(
        description='Enhance cimage html output with preview script'
    )

    parser.add_argument(
        '-b', '--backup', 
        help='Backup file before adding script',
        action='store_true'
    )

    parser.add_argument(
        '-q', '--quick',
        help='Do simple string substitution without parsing html',
        action='store_true'
    )

    parser.add_argument(
        '-r', '--remove',
        help='Removes script from files',
        action='store_true'
    )

    parser.add_argument(
        '-R', '--removeall',
        help='Removes ALL scripts from files',
        action='store_true'
    )

    parser.add_argument('args', nargs=REMAINDER)

    return parser.parse_args()

if __name__ == '__main__':
    main()
